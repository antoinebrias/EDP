function mdlstruct = fit_gp(mdlstruct)
% FIT_GP   Gaussian Process fitting.
%    MDLSTRUCT = FIT_GP(MDLSTRUCT) normalizes the data in MDLSTRUCT and
%    optimize hyperparameters of a Gaussian Process regression on these
%    data, according the constraints available in the MDLSTRUCT.GP
%    structure:
%
%    GP.IS_LOG indicates if we consider fitting log(X_{t+1}/X_t) instead
%    of X_{t+1}
%
%    GP.COND0 indicates if a constraint on [0,0,...,0] is applied.
%
%    GP.COND_J_0 indicates if a constraint on 0 is applied on the j-th
%    variable.
%
%    GP.COND_J_1 indicates if a constraint on 1 is applied on the j-th
%    variable. !!! Don't apply more than one constraint for now.
%
%   As output, you will find:
%    GP.GP_MODEL is a handle of the posterior GP.
%    GP.IS_ERR is the in-sample error of the posterior GP.
%    GP.OOS_ERR is the out-of-sample error of the posterior GP.
%    GP.LENGTHSCALES is the lenghscale parameters for each variables.
%    GP.NLL is the negative-likelihood.
%    GP.PFIT is the hyperparameter vector of the Gaussian Process
%
%   See also POST_GP, FIT_VALUE_FUNCTION.



gp = mdlstruct.gp;

is_log = gp.is_log;
cond0 = gp.cond_0;
Condj0=gp.cond_j_0;
Condj1 =gp.cond_j_1;

 mdlstruct.gp.is_value_function=0;

n_dim = mdlstruct.n_dim;
tau = mdlstruct.tau;

n_lags = mdlstruct.n_lags; % number of lags for each covariables

% Index of first column of each species (used with TD dim)
ind_available_var=find(n_lags>0);
ind_current_var=unique(cumsum(n_lags)-n_lags+ones(1,length(n_lags)));
ind_current_var(ind_current_var>sum(n_lags))=[];

Data  =mdlstruct.data(:,ind_available_var);
Control = mdlstruct.control_data(:,ind_available_var);


switch mdlstruct.control_type
    case 'rate'
        if iscell(mdlstruct.control_param)
            
        else
            for i = 1:length(ind_available_var) % !!!14/04/2021 change n_dim
                %                 Input(:,i) = Data((n_lags-1)*tau+1:end-tau,i).*(1-Control((n_lags-1)*tau+1:end-tau,i)*mdlstruct.control_param(i));
                %                 Output(:,i)= Data(n_lags*tau+1:end,i);
                
%                 Input(:,i) = Data((n_lags(ind_available_var(i))-1)*tau:end-tau,ind_available_var(i)).*(1-Control((n_lags(ind_available_var(i))-1)*tau:end-tau,ind_available_var(i))*mdlstruct.control_param(ind_available_var(i)));
%               
%                 Output(:,i)= Data(n_lags(ind_available_var(i))*tau:end,ind_available_var(i));
%                 
                Input(:,i) = Data(1:end-1,ind_available_var(i)).*(1-Control(1:end-1,ind_available_var(i))*mdlstruct.control_param(ind_available_var(i)));
                Output(:,i)= Data(2:end,ind_available_var(i));
            end
        end
        
    case 'single'
        
    case 'global'
        
end


% data normalization
mIn=0.*min(Input);
sdIn=1+0.*(max(Input)-1*min(Input));


% %% Change 25/11/2021
% s_in = Input;
% y = Output;



s_in = Input;
y = Output;


%% Lag matrix construction
x = [];
for i=1:length(ind_available_var) % !!!14/04/2021 change n_dim
    %lag matrix
    xtmp = hankel_matrix((Input(:,i)-mIn(i))/(sdIn(i)+eps),[0:tau:tau*(n_lags(ind_available_var(i))-1)]);
    x=[x xtmp(:,1:end)];
end
% %remove extra rows
if n_lags==1
y(1, :) = [];
s_in(1, :) = [];
x(1, :) = [];
end
y(end, :) = [];
s_in(end, :) = [];
x(end, :) = [];
% if n_lags==2
% y(end, :) = [];
% s_in(end, :) = [];
% x(end, :) = [];
% end


%remove rows with nan
y(any(isnan(x), 2), :) = [];
s_in(any(isnan(x), 2), :) = [];
x(any(isnan(x), 2), :) = [];


if is_log
    z = log(y./(s_in+eps)+eps);%log((Data(nlags*tau+1:end,:)+eps)./(Data((nlags-1)*tau+1:end-tau,:)+eps)+eps); % z = log(y/s)=f(s).
else
    z =y;
end

%%%% Change 25/11/2021
if  ~mdlstruct.gp.is_value_function
    z = log(y);%!!!!!!!!!!!!!!!!!!!!
x = log(x+eps);
y = log(y+eps);
s_in = log(s_in+eps);

end




%
% mOut=min(z);%min(z); %%%%%% CHANGE 18/11/21
% sdOut=std(z);

mOut=mean(z);%min(z); %%%%%% CHANGE 25/11/21
sdOut=std(z);


zOut=(z-mOut)./(sdOut+eps);

%




%% GP hyperparameters optim;
LenScale=.05*ones(1,sum(n_lags));
ve=.000000001;
gp_tau=10;
pars=[log(LenScale) log(ve/(1-ve)) log(gp_tau)]';

%     cond0 = 0;

%




for i=1:length(ind_available_var) % !!!14/04/2021 change n_dim
    
    lpost=@(p) GP4DP(p,x,zOut(:,i),[],[],cond0,Condj0,Condj1,0);
    [pfit,nll]=fmingrad_Rprop(lpost,pars);
    phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
    Pfit(:,i) = pfit;
    GP{i}= @(X)GP4DP(pfit,x,zOut(:,i),X,[],cond0,Condj0,Condj1,0); %handle function x does the long function
    Nll(i) = nll;
    %     end
    
    [dummy1,dummy2,out]=GP{i}(x);
    abserr=abs(out.pred'-zOut(:,i)');
    errIS(i)=mean(abserr.^2)./(var(zOut(:,i))+eps); % err shows err values
    
    
    for j=1:length(z(:,i))
        xnoj=x;xnoj(j,:)=[];
        ynoj=zOut(:,i);ynoj(j)=[];
        xj=x(j,:);yj=zOut(j,i);
        %     [dummy1,dummy2,out]=GP4DP_Embed2b_fit(pfit,xnoj,ynoj,[],xj,nlags*2, 0,condpars);
        [dummy1,dummy2,out]=GP4DP(pfit,xnoj,ynoj,xj,[],cond0,Condj0,Condj1,0);
        yjpred(j)=out.pred;
    end
    
    abserr=abs(yjpred-zOut(:,i)');
    errOS(i)=mean(abserr.^2)./(var(zOut(:,i))+eps); % err shows err values
    
    
end

mdlstruct.ind_available_var =ind_available_var;
mdlstruct.ind_current_var=ind_current_var;

mdlstruct.gp.ind_available_var=ind_available_var;
mdlstruct.gp.is_err=errIS;
mdlstruct.gp.oos_err=errOS;
mdlstruct.gp.gp_handle=GP;
mdlstruct.gp.n_dim=n_dim;
mdlstruct.gp.tau=tau;
mdlstruct.gp.n_lags=n_lags;
mdlstruct.gp.lengthscales=phi;
mdlstruct.gp.pfit=Pfit;
mdlstruct.gp.nll=Nll;
mdlstruct.gp.mIn=mIn;
mdlstruct.gp.sdIn=sdIn;
mdlstruct.gp.mOut=mOut;
mdlstruct.gp.sdOut=sdOut;



% mdlstruct.gp_model = @(x,u,is_det,is_hist)post_gp(x,u,is_det,is_hist,mdlstruct.gp);
mdlstruct.gp_model = @(x,u,is_det)post_gp(x,u,is_det,mdlstruct.gp);

end


%
%
% %
% for i=1:nvar
%     %     if cond0
%     %     lpost=@(p) GP4DP(p,s,z(:,i),[],[],0,[i 0],Condj1,0);
%     %     [pfit,nll]=fmingrad_Rprop(lpost,pars);
%     %     phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
%     %     Pfit(:,i) = pfit;
%     %     GPi{i}= @(X)GP4DP(pfit,s,z(:,i),X,[],0,[i 0],Condj1,0); %handle function x does the long function
%     %     Nll(i) = nll;
%     %     else
%
%     lpost=@(p) GP4DP(p,s,z(:,i),[],[],cond0,Condj0,Condj1,0);
%     [pfit,nll]=fmingrad_Rprop(lpost,pars);
%     phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
%     Pfit(:,i) = pfit;
%     GPi{i}= @(X)GP4DP(pfit,s,z(:,i),X,[],cond0,Condj0,Condj1,0); %handle function x does the long function
%     Nll(i) = nll;
%     %     end
%
%     [dummy1,dummy2,out]=GPi{i}(s);
%     abserr=abs(out.pred'-z(:,i)');
%     errIS(i)=mean(abserr.^2)./(var(z(:,i))+eps); % err shows err values
% end
% GP.GPi = GPi;
% GP.mOut = mOut;    GP.mIn = mIn;
% GP.sdOut = sdOut;  GP.sdIn = sdIn;
% GP.isLog = isLog;
% GP.nSpecies = nvar;
% GP.Pfit=Pfit;
% GP.Nll=Nll;
%
% if isDebug
%     %% out of sample loop
%     yjpred=[];
%     for i=1:nvar
%         for j=1:length(z(:,i))
%             xnoj=s;xnoj(j,:)=[];
%             ynoj=z(:,i);ynoj(j)=[];
%             xj=s(j,:);yj=z(j);
%             %     [dummy1,dummy2,out]=GP4DP_Embed2b_fit(pfit,xnoj,ynoj,[],xj,nlags*2, 0,condpars);
%             [dummy1,dummy2,out]=GP4DP(pfit,xnoj,ynoj,xj,[],cond0,Condj0,Condj1,0);
%             yjpred(j)=out.pred;
%         end
%
%         abserr=abs(yjpred-z(:,i)');
%         errOS(i)=mean(abserr.^2)./(var(z(:,i))+eps); % err shows err values
%         outOfSamplePred = yjpred;
%     end
%     errIS
%     errOS
% %
% %
% % function [GP]=fitGP(input,output,nSpecies, isLog, cond0)
% % %% data treatment
% % Data  = input;
% %
% % % data normalization
% % mIn=min(Data)*1;
% % sdIn=max(Data)-1*min(Data);
% %
% %
% % sIn=(Data-mIn)./(sdIn+eps);
% % sData = sIn;
% %
% % maxOut = max(output);
% %
% % if isLog
% %     z = log(output./(input+eps)+eps); % z = log(y/s)=f(s).
% % else
% %     z =output;
% % end
% %
% % mOut=1*min(z);
% % sdOut=max(z)-1*min(z);
% %
% % mOut=min(z);
% % sdOut=std(z);
% % % if isLog
% % %     mOut=0*min(z)-0*10;%min(z)*0;;
% % %     sdOut=max(z)-min(z)+0*10;
% % %
% % % %         mOut=mean(z);%min(z)*0;;
% % % %     sdOut=std(z);
% % %
% % % end
% % zOut=(z-mOut)./(sdOut+eps);
% % zData = zOut;
% %
% % % mOut = mOut.*0;
% %
% % nlags = 1; % number of lags for each covariables
% % tau = 1;  % tau
% %
% %
% % z = zData;
% % s = sData;
% %
% % %% GP hyperparameters optim;
% % condpars=[0 0];
% % LenScale=.05*ones(1,nlags*nSpecies);
% % ve=.00000001;
% % tau=10;
% % pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';
% %
% % %     cond0 = 0;
% % Condj0=[0 0];Condj1 =[0 0];
% % %
% % %
% % for i=1:nSpecies
% %     %     if cond0
% %     %     lpost=@(p) GP4DP(p,s,z(:,i),[],[],0,[i 0],Condj1,0);
% %     %     [pfit,nll]=fmingrad_Rprop(lpost,pars);
% %     %     phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
% %     %     Pfit(:,i) = pfit;
% %     %     GPi{i}= @(X)GP4DP(pfit,s,z(:,i),X,[],0,[i 0],Condj1,0); %handle function x does the long function
% %     %     Nll(i) = nll;
% %     %     else
% %
% %     lpost=@(p) GP4DP(p,s,z(:,i),[],[],cond0,Condj0,Condj1,0);
% %     [pfit,nll]=fmingrad_Rprop(lpost,pars);
% %     phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
% %     Pfit(:,i) = pfit;
% %     GPi{i}= @(X)GP4DP(pfit,s,z(:,i),X,[],cond0,Condj0,Condj1,0); %handle function x does the long function
% %     Nll(i) = nll;
% %     %     end
% % end
% % GP.GPi = GPi;
% % GP.mOut = mOut;    GP.mIn = mIn;
% % GP.sdOut = sdOut;  GP.sdIn = sdIn;
% % GP.isLog = isLog;
% % GP.nSpecies = nSpecies;
% % GP.Pfit=Pfit;
% % GP.Nll=Nll;
% % GP.maxOut = maxOut;
% % end
% %
% %
% %
% % function [mu,v]=postMdl(Mdl,s,noNoise,isValueFct)
% % if iscell(Mdl)
% %     for i=1:length(Mdl)
% %         mu(:,i)=predict(Mdl{i},s(:,i));
% %     end
% %     v = mu*0;
% % else
% %
% %     % single species policy value %% for init of td algo
% %     if  isfield(Mdl,'vpred')
% %         % value of 1sp policy (unweighted)
% %         vrespreys = [predict(Mdl.vpreys{1},sum(s(:,[1 2 3]),2))];
% %         vrespred = [predict(Mdl.vpred{1},s(:,4))];
% %         mu = [vrespreys./3 vrespreys./3 vrespreys./3 vrespred zeros(size(s,1),1)];
% %         v = mu*0;
% %     else
% %
% %         snorm = (s-Mdl.mIn)./(Mdl.sdIn);
% %
% %         for i=1:Mdl.nSpecies
% %             [~,~,out]=Mdl.GPi{i}(snorm);
% %             muTmp(:,i)=out.pred*Mdl.sdOut(i)+Mdl.mOut(i);
% %             if Mdl.isLog==1
% %                 mu(:,i) = exp( muTmp(:,i)).*s(:,i);
% %             else
% %                 mu(:,i) =  muTmp(:,i);
% %             end
% %             v(:,i) = out.var*Mdl.sdOut(i)^2;
% %
% %
% %             if  nargin> 3 &  isValueFct
% %                 % avoid issue near 0
% %                 mu(s(:,i)<repmat(0.5*Mdl.mIn(:,i),size(mu(:,i),1),1),i)=0;
% %                 %% avoid issue when GP post mean is higher than data max
% %                 mu(mu(:,i)>repmat(1.1*Mdl.maxOut(:,i),size(mu(:,i),1),1),i)=1.1*Mdl.maxOut(:,i);
% %
% %
% %                 if isfield(Mdl,'isLinearGP') & Mdl.isLinearGP==1
% %                     tmp = exp(log(s+eps)*Mdl.linearCoef);
% %                     mu(:,i) = mu(:,i) +  tmp(:,i);
% %                 end
% %
% %
% %
% %
% %             end
% %
% %
% %         end
% %     end
% %
% %     if nargin> 2 & ~noNoise
% %         if Mdl.isLog==1
% %             mu = s.*exp(normrnd(muTmp,sqrt(v)));
% %         else
% %             mu = normrnd(muTmp,sqrt(v));
% %         end
% %     end
% % end
% % mu(mu<0)=0;
% %
% % % if nargin> 3 & isValueFct
% %
% % end
