function value_function = fit_value_function(value_function,input,output,mdlstruct)
% FIT_VALUE_FUNCTION   Gaussian Process regression for value function.
%    VALUE_FUNCTION = FIT_VALUE_FUNCTION(VALUE_FUNCTION,INPUT,OUTPUT,MDLSTRUCT) fits
%    a Gaussan Process on the states in INPUT with their value in OUTPUT.
%    
%   As output, you will find:
%    VALUE_FUNCTION.GP_MODEL is a handle of the posterior GP.
%    VALUE_FUNCTION.IS_ERR is the in-sample error of the posterior GP.
%    VALUE_FUNCTION.OOS_ERR is the out-of-sample error of the posterior GP.
%    VALUE_FUNCTION.LENGTHSCALES is the lenghscale parameters for each variables.
%    VALUE_FUNCTION.NLL is the negative-likelihood.
%    VALUE_FUNCTION.PFIT is the hyperparameter vector of the Gaussian
%    Process
%
%   See also FIT_GP.


if isfield(value_function,'gp')
gp = value_function.gp;
else
   gp = []; 
end

if isfield(gp,'is_log')
    % TD lambda parameter
    is_log = gp.is_log;
else
    is_log=0;
end

if isfield(gp,'cond_0')
    % TD lambda parameter
    cond_0 = gp.cond_0;
else
    cond_0=0;
end

if isfield(gp,'cond_j_0')
    % TD lambda parameter
    cond_j_0 = gp.cond_j_0;
else
    cond_j_0=[0 0];
end


if isfield(gp,'cond_j_1')
    % TD lambda parameter
    cond_j_1 = gp.cond_j_1;
else
    cond_j_1=[0 0];
end


n_dim_input = size(input,2);
n_dim_output = size(output,2);
tau = 1;

nlags = 1; % number of lags for each covariables


n_dim = mdlstruct.n_dim;
tau = mdlstruct.tau;

n_lags = mdlstruct.n_lags; % number of lags for each covariables

% Index of first column of each species (used with TD dim)
ind_available_var=find(n_lags>0);
ind_current_var=unique(cumsum(n_lags)-n_lags+ones(1,length(n_lags)));
ind_current_var(ind_current_var>sum(n_lags))=[];

% data normalization
mIn=min(input)*1;
sdIn=max(input)-1*min(input);


%%
s_in = input;
y = output;
x = (input - mIn)./(sdIn+eps);

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


%% GP hyperparameters optim;
LenScale=.05*ones(1,nlags*n_dim_input);
ve=.00000001;
gp_tau=10;
pars=[log(LenScale) log(ve/(1-ve)) log(gp_tau)]';

%     cond0 = 0;
for i=1:n_dim_output
    
    lpost=@(p) GP4DP(p,x,zOut(:,i),[],[],cond_0,cond_j_0,cond_j_1,0);
    [pfit,nll]=fmingrad_Rprop(lpost,pars);
    phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
    Pfit(:,i) = pfit;
    GP{i}= @(X)GP4DP(pfit,x,zOut(:,i),X,[],cond_0,cond_j_0,cond_j_1,0); %handle function x does the long function
    Nll(i) = nll;
    %     end
end

value_function.gp.is_log=is_log;
value_function.gp.cond_0=cond_0;
value_function.gp.cond_j_0=cond_j_0;
value_function.gp.cond_j_1=cond_j_1;

value_function.gp.ind_available_var=ind_available_var;
value_function.gp.gp_handle=GP;
value_function.gp.n_dim=n_dim_output;
value_function.gp.lengthscales=phi;
value_function.gp.pfit=Pfit;
value_function.gp.nll=Nll;
value_function.gp.mIn=mIn;
value_function.gp.sdIn=sdIn;
value_function.gp.mOut=mOut;
value_function.gp.sdOut=sdOut;
value_function.gp_model = @(x,is_det)post_gp(x,0,is_det,value_function.gp);

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
