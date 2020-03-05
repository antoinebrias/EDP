%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Approximate dynamic programming algorithm
% Using temporal difference learning
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dpstruct] = td_learning(mdlstruct,optstruct,dpstruct,name_model)



%% optimal control problem parameters
% Discount factor
discount_factor=optstruct.discount_factor;

% Objective weights
weights=optstruct.weights;

% Reward function
reward=optstruct.reward;

%% Temporal difference learning parameters
% number of support states
n_support_states=dpstruct.n_support_states;


% TD lambda parameter
lambda=dpstruct.lambda;

% outer loop number of iterations
n_out_max=dpstruct.n_out_max;

% inner loop number of iterations
n_in_max=dpstruct.n_in_max;

% debug mode with additionnal display during the computation
debug=dpstruct.debug;

if ~isfield(dpstruct,'value_function_init')
        value_function=[];
else
    value_function = dpstruct.value_function_init;
end


 
if ~isfield(dpstruct,'support_states')
    [dpstruct] = generate_support_states(mdlstruct,dpstruct);
else
       support_states = dpstruct.support_states;
end


switch name_model
    case 'gp'
        model=mdlstruct.gp_model;
    case 'mdl'
        model= mdlstruct.model;
    otherwise
        
end


% control grid
U=linspace(0,1,10);
[~,uWeightsInd]=find(mdlstruct.control_param>0);
switch mdlstruct.control_type
    case 'rate'
        if iscell(mdlstruct.control_param)
            % to do
        else
             
            [gridTmpU{1:nnz(mdlstruct.control_param)}] = deal(U);
            [gridTmpU{1:nnz(mdlstruct.control_param)}] = ndgrid(gridTmpU{:});
            control_grid=[];
            for i=1:length(mdlstruct.control_param)
                [isMember,ii]=ismember(i,uWeightsInd);
                if isMember
                    control_grid=[control_grid gridTmpU{ii}(:)];
                else
                    control_grid=[control_grid zeros(nU^nnz(uWeights),1)];
                    % columns of zeros for uncontroled species
                end
            end
            
            
            
            
        end
    case 'single'
        [gridTmpU{1}] = deal(U);
        [gridTmpU{1}] = ndgrid(gridTmpU{:});
        control_grid=[];
        for i=1:length(uWeights)
            [isMember,ii]=ismember(i,uWeightsInd);
            if isMember
                control_grid=[control_grid gridTmpU{1}(:).*mdlstruct.control_param(i)];
            else
                control_grid=[control_grid zeros(nU,1)];
            end
        end
        control_grid= unique(control_grid,'rows');
    case 'global'
        % to do
    otherwise
        
end



n_dim = mdlstruct.n_dim;
ssInit = dpstruct.support_states;


for n_out=1:n_out_max
    ss = ssInit;
    for n_in=1:n_in_max
        currentX =ss;
        
        %% One-step ahead DP problem
        % solve the one-step ahead DP equation by looking for the control
        % maximizing fun_optim_td
        fOptTD = @(control) fun_optim_td (support_states,control,value_function,model,optstruct,dpstruct);
        
        ftd = [];
        
        [ftd,fitmp]=fOptTD(control_grid);
        
        
        [fTD,indCurrentU]=min(reshape(ftd,length(control_grid),length(currentX)));
        
        % optimal control
        currentU=control_grid(indCurrentU,:);
        
        fitmp2=reshape(fitmp,length(control_grid),length(currentX),n_dim);
        for i=1:length(currentX)
            fI(:,i) =squeeze(fitmp2(indCurrentU(i),i,:));
        end
        
        
        %% Computing temporal difference
        
        if n_out==1 & isempty(value_function)
            [vTD(:,n_in),vI(:,:,n_in)]=  reward(ss,ss.*0);
        else
            [muV,vV]=value_function.gp_model(ss,1);
            vI(:,:,n_in)=muV;
            vI(vI<0)=0;
            
        end
        
        % Temporal difference
        Di=-fI' -vI(:,:,n_in);
        
        % Several learning rate are possible
        if n_in==1
            alphaInner = 1;
        else
            alphaInner=5/(5+n_in-2); % harmonic
        end
        totDiscount = (optstruct.discount_factor*lambda).^(n_in-(1:n_in));
        
        %% Step 4 - Update and store the value function + store each species part of the value (vI)
        vI(:,:,1:n_in)=vI(:,:,1:n_in)+alphaInner*Di.*reshape(totDiscount,1,[],length(totDiscount));
        
        vI(vI<0)=0;
        
        %% Simulating trajectories
        [mu,~]=model(ss,currentU,0 );
        nextX = mu;
        nextX(nextX<0 )=0;
        
        ssNew = nextX;
        ss=ssNew;
        
        %         if n_in ==1
        %             Usave =currentU;
        %         end
    end
    
    %%   Value function approximation
    
    if n_out==1
        
        value_function = fit_value_function(value_function,ssInit,vI(:,:,1));
        alphaOuter=1;
        
        
    else
        %We update the GP according the new values of the support states
        alphaOuter=1;%5/(5+nOut-2); % harmonic
        [muV,vV]=value_function.gp_model(ssInit,1);
        vPrev = muV; vIPrev=muV;
        vCurrent = ((1-alphaOuter)*vIPrev+alphaOuter*vI(:,:,1));
        
        
        value_function = fit_value_function(value_function,ssInit,vCurrent);
        
    end
    
    % %
    %                         figure(1123);
    %                         subplot(2,2,1)
    %                         plot(ssInit(:,1),vI(:,:,1),'.','Linewidth',2)
    %                         hold on
    %                         sss = repmat(mean(ssInit),50,1);
    %                         xx = linspace(0,max(ssInit(:,1)),50);sss(:,1)=xx';
    %                         [mu,v]= postGP(GPV,sss,[],1);hold on; plot(xx,mu,'r')
    %                         hold off
    %                         subplot(2,2,2)
    %                         plot(ssInit(:,1),Usave(:,1),'.','Linewidth',2)
    %                         %      subplot(2,2,3)
    %                         %     plot(paramTD.ssInit.currentX{1},vI(:,:,1).*paramReward,'.','Linewidth',2)
    %                         subplot(2,2,[3 4])
    %                         plot(aaa);
    %                         pause(0.01)
    %                             nOut
    
end

dpstruct.model = model;
dpstruct.control_grid=control_grid;
dpstruct.value_function=value_function;
end







