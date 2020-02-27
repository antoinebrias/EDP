%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gaussian Process regression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mdlstruct = fit_gp(mdlstruct)
gp = mdlstruct.gp;
input = gp.mdl_type
output = gp.mdl_type

isLog = gp.is_log
cond0 = gp.cond_0;
Condj0=gp.cond_j_0;
Condj1 =gp.cond_j_1;
ndim = mdlstruct.n_dim;
tau = mdlstruct.tau;

nlags = mdlstruct.n_lags_max; % number of lags for each covariables



nTD =nlags*ndim;




Data  =mdlstruct.data;

% data normalization
mdata=min(Data).*0;
sddata=max(Data)-0*min(Data);
dataNorm=(Data-mdata)./(sddata+eps);

xData = dataNorm;
yData = dataNorm;


tau = 1;  % tau

%% Lag matrix construction
x = [];
y = yData;
for i=1:2
    %lag matrix
    xtmp = lagmatrix(xData(:,i),[0:tau:tau*nlags]);
    x=[x xtmp(:,2:end)];
end
%remove rows with nan
y(any(isnan(x), 2), :) = [];
x(any(isnan(x), 2), :) = [];


%% GP hyperparameters optim;
LenScale=.05*ones(1,nlags*2);
ve=.000000001;
Tau=10;
pars=[log(LenScale) log(ve/(1-ve)) log(Tau)]';




lpost=@(p) GP4DP(p,x,y,[],[],cond0,Condj0,Condj1,0);
[pfit,nll]=fmingrad_Rprop(lpost,pars);
phi = exp(pfit(1:end-2))'; % lengthscale parameters


% GP= @(X)GP4DP_Embed2b_fit(pfit,x,y,[],X,nlags*2, 0,condpars); %handle function x does the long function
GP= @(X)GP4DP(pfit,x,y,X,[],cond0,Condj0,Condj1,0); %handle function x does the long function
[neglpost,neglgrad,out]=GP(x); %out.pred and out.var give the posterior mean and var






mdlstruct.gp.lengthscales=phi;
mdlstruct.gp.pfit=pfit;
mdlstruct.gp.nll=nll;
mdlstruct.model = @(x,u,);



end



function [GP]=fitGP(input,output,nSpecies, isLog, cond0)
%% data treatment
Data  = input;

% data normalization
mIn=min(Data)*1;
sdIn=max(Data)-1*min(Data);


sIn=(Data-mIn)./(sdIn+eps);
sData = sIn;

maxOut = max(output);

if isLog
    z = log(output./(input+eps)+eps); % z = log(y/s)=f(s).
else
    z =output;
end

mOut=1*min(z);
sdOut=max(z)-1*min(z);

mOut=min(z);
sdOut=std(z);
% if isLog
%     mOut=0*min(z)-0*10;%min(z)*0;;
%     sdOut=max(z)-min(z)+0*10;
%
% %         mOut=mean(z);%min(z)*0;;
% %     sdOut=std(z);
%
% end
zOut=(z-mOut)./(sdOut+eps);
zData = zOut;

% mOut = mOut.*0;

nlags = 1; % number of lags for each covariables
tau = 1;  % tau


z = zData;
s = sData;

%% GP hyperparameters optim;
condpars=[0 0];
LenScale=.05*ones(1,nlags*nSpecies);
ve=.00000001;
tau=10;
pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';

%     cond0 = 0;
Condj0=[0 0];Condj1 =[0 0];
%
%
for i=1:nSpecies
    %     if cond0
    %     lpost=@(p) GP4DP(p,s,z(:,i),[],[],0,[i 0],Condj1,0);
    %     [pfit,nll]=fmingrad_Rprop(lpost,pars);
    %     phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
    %     Pfit(:,i) = pfit;
    %     GPi{i}= @(X)GP4DP(pfit,s,z(:,i),X,[],0,[i 0],Condj1,0); %handle function x does the long function
    %     Nll(i) = nll;
    %     else
    
    lpost=@(p) GP4DP(p,s,z(:,i),[],[],cond0,Condj0,Condj1,0);
    [pfit,nll]=fmingrad_Rprop(lpost,pars);
    phi(i,:) = exp(pfit(1:end-2))'; % lengthscale parameters
    Pfit(:,i) = pfit;
    GPi{i}= @(X)GP4DP(pfit,s,z(:,i),X,[],cond0,Condj0,Condj1,0); %handle function x does the long function
    Nll(i) = nll;
    %     end
end
GP.GPi = GPi;
GP.mOut = mOut;    GP.mIn = mIn;
GP.sdOut = sdOut;  GP.sdIn = sdIn;
GP.isLog = isLog;
GP.nSpecies = nSpecies;
GP.Pfit=Pfit;
GP.Nll=Nll;
GP.maxOut = maxOut;
end
