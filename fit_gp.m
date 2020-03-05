%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gaussian Process regression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mdlstruct = fit_gp(mdlstruct)
gp = mdlstruct.gp;


is_log = gp.is_log;
cond0 = gp.cond_0;
Condj0=gp.cond_j_0;
Condj1 =gp.cond_j_1;



n_dim = mdlstruct.n_dim;
tau = mdlstruct.tau;

nlags = mdlstruct.n_lags_max; % number of lags for each covariables

nTD =nlags*n_dim;

Data  =mdlstruct.data;
Control = mdlstruct.control_data;


switch mdlstruct.control_type
    case 'rate'
        if iscell(mdlstruct.control_param)
            
        else
            for i = 1:length(mdlstruct.control_param)
                Input(:,i) = Data((nlags-1)*tau+1:end-tau,i).*(1-Control((nlags-1)*tau+1:end-tau,i)*mdlstruct.control_param(i));
                Output(:,i)= Data(nlags*tau+1:end,i);
                %                 Input(:,i) = Data(:,i).*(1-Control(:,i));
                %                 Output(:,i)= Data(:,i);
            end
        end
        
    case 'single'
        
    case 'global'
        
end


% data normalization
mIn=min(Input)*1;
sdIn=max(Input)-1*min(Input);


%%
s_in = Input;
y = Output;
%% Lag matrix construction
x = [];
for i=1:n_dim
    %lag matrix
    xtmp = lag_matrix((Input(:,i)-mIn(i))/(sdIn(i)+eps),[0:tau:tau*(nlags-1)]);
    x=[x xtmp(:,1:end)];
end
%remove rows with nan
y(any(isnan(x), 2), :) = [];
s_in(any(isnan(x), 2), :) = [];
x(any(isnan(x), 2), :) = [];


if is_log
    z = log(y./(s_in+eps)+eps);%log((Data(nlags*tau+1:end,:)+eps)./(Data((nlags-1)*tau+1:end-tau,:)+eps)+eps); % z = log(y/s)=f(s).
else
    z =y;
end

%
mOut=min(z);
sdOut=std(z);




zOut=(z-mOut)./(sdOut+eps);
zData = zOut;
%




%% GP hyperparameters optim;
LenScale=.05*ones(1,nlags*n_dim);
ve=.00000001;
gp_tau=10;
pars=[log(LenScale) log(ve/(1-ve)) log(gp_tau)]';

%     cond0 = 0;

%
%
for i=1:n_dim
    
    lpost=@(p) GP4DP(p,x,zOut(:,i),[],[],cond0,Condj0,Condj1,0);
    [pfit,nll]=fmingrad_Rprop(lpost,pars);
    phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
    Pfit(:,i) = pfit;
    GP{i}= @(X)GP4DP(pfit,x,zOut(:,i),X,[],cond0,Condj0,Condj1,0); %handle function x does the long function
    Nll(i) = nll;
    %     end
end

mdlstruct.gp.gp_handle=GP;
mdlstruct.gp.n_dim=n_dim;
mdlstruct.gp.lengthscales=phi;
mdlstruct.gp.pfit=Pfit;
mdlstruct.gp.nll=Nll;
mdlstruct.gp.mIn=mIn;
mdlstruct.gp.sdIn=sdIn;
mdlstruct.gp.mOut=mOut;
mdlstruct.gp.sdOut=sdOut;
mdlstruct.gp_model = @(x,u,is_det)post_gp(x,u,is_det,mdlstruct.gp);

end

%
%
% function [GP]=fitGP(input,output,nSpecies, isLog, cond0)
% %% data treatment
% Data  = input;
%
% % data normalization
% mIn=min(Data)*1;
% sdIn=max(Data)-1*min(Data);
%
%
% sIn=(Data-mIn)./(sdIn+eps);
% sData = sIn;
%
% maxOut = max(output);
%
% if isLog
%     z = log(output./(input+eps)+eps); % z = log(y/s)=f(s).
% else
%     z =output;
% end
%
% mOut=1*min(z);
% sdOut=max(z)-1*min(z);
%
% mOut=min(z);
% sdOut=std(z);
% % if isLog
% %     mOut=0*min(z)-0*10;%min(z)*0;;
% %     sdOut=max(z)-min(z)+0*10;
% %
% % %         mOut=mean(z);%min(z)*0;;
% % %     sdOut=std(z);
% %
% % end
% zOut=(z-mOut)./(sdOut+eps);
% zData = zOut;
%
% % mOut = mOut.*0;
%
% nlags = 1; % number of lags for each covariables
% tau = 1;  % tau
%
%
% z = zData;
% s = sData;
%
% %% GP hyperparameters optim;
% condpars=[0 0];
% LenScale=.05*ones(1,nlags*nSpecies);
% ve=.00000001;
% tau=10;
% pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';
%
% %     cond0 = 0;
% Condj0=[0 0];Condj1 =[0 0];
% %
% %
% for i=1:nSpecies
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
% end
% GP.GPi = GPi;
% GP.mOut = mOut;    GP.mIn = mIn;
% GP.sdOut = sdOut;  GP.sdIn = sdIn;
% GP.isLog = isLog;
% GP.nSpecies = nSpecies;
% GP.Pfit=Pfit;
% GP.Nll=Nll;
% GP.maxOut = maxOut;
% end
%
%
%
% function [mu,v]=postMdl(Mdl,s,noNoise,isValueFct)
% if iscell(Mdl)
%     for i=1:length(Mdl)
%         mu(:,i)=predict(Mdl{i},s(:,i));
%     end
%     v = mu*0;
% else
%
%     % single species policy value %% for init of td algo
%     if  isfield(Mdl,'vpred')
%         % value of 1sp policy (unweighted)
%         vrespreys = [predict(Mdl.vpreys{1},sum(s(:,[1 2 3]),2))];
%         vrespred = [predict(Mdl.vpred{1},s(:,4))];
%         mu = [vrespreys./3 vrespreys./3 vrespreys./3 vrespred zeros(size(s,1),1)];
%         v = mu*0;
%     else
%
%         snorm = (s-Mdl.mIn)./(Mdl.sdIn);
%
%         for i=1:Mdl.nSpecies
%             [~,~,out]=Mdl.GPi{i}(snorm);
%             muTmp(:,i)=out.pred*Mdl.sdOut(i)+Mdl.mOut(i);
%             if Mdl.isLog==1
%                 mu(:,i) = exp( muTmp(:,i)).*s(:,i);
%             else
%                 mu(:,i) =  muTmp(:,i);
%             end
%             v(:,i) = out.var*Mdl.sdOut(i)^2;
%
%
%             if  nargin> 3 &  isValueFct
%                 % avoid issue near 0
%                 mu(s(:,i)<repmat(0.5*Mdl.mIn(:,i),size(mu(:,i),1),1),i)=0;
%                 %% avoid issue when GP post mean is higher than data max
%                 mu(mu(:,i)>repmat(1.1*Mdl.maxOut(:,i),size(mu(:,i),1),1),i)=1.1*Mdl.maxOut(:,i);
%
%
%                 if isfield(Mdl,'isLinearGP') & Mdl.isLinearGP==1
%                     tmp = exp(log(s+eps)*Mdl.linearCoef);
%                     mu(:,i) = mu(:,i) +  tmp(:,i);
%                 end
%
%
%
%
%             end
%
%
%         end
%     end
%
%     if nargin> 2 & ~noNoise
%         if Mdl.isLog==1
%             mu = s.*exp(normrnd(muTmp,sqrt(v)));
%         else
%             mu = normrnd(muTmp,sqrt(v));
%         end
%     end
% end
% mu(mu<0)=0;
%
% % if nargin> 3 & isValueFct
%
% end
