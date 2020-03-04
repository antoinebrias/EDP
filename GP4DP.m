function[neglpost,neglgrad,out,dataGP]=GP4DP(pars,X,Y,predgrid,dataGP,cond0,Condj0,Condj1,noVe)
%This version assumes X is scaled [0,1] and Y is centered
%covariance function is
%exponential (dexp=1), squared-exponential (dexp=2), or conditional on f(0,0)=0 (cond0=1);
%Condj0 : condition when f(x|xj=0)=Fcondj0
%Condj1 : condition when f(x|xj=1)=Fcondj1
if nargin <7
    Condj0(1) = 0;
    condj0=0;
    Fcondj0=0;
           Condj1(1) = 0;
        condj1=0;
        Fcondj1=0;
else
    condj0 = Condj0(1)>0;
    Fcondj0 =  Condj0(2);
    if nargin <8
        Condj1(1) = 0;
        condj1=0;
        Fcondj1=0;
    else
        condj1 = Condj1(1)>0;
        Fcondj1 =  Condj1(2);
    end
end
h=1;
T=length(X);
dexp=2;
%dexp=1;

d=size(X,2);
npars=d+2;

vec=@(x)x(:);
%pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';




%transform parameters from real line to constrained space
vemin=0.001;taumin=.001;
phi=exp(pars(1:d));%phi=0.1*ones(d,1);
ve=(1-vemin)*exp(pars(d+1))/(1+exp(pars(d+1)))+vemin;
tau=(1-taumin)*exp(pars(d+2))/(1+exp(pars(d+2)))+taumin;
dpars=[phi;ve*(1-vemin-ve)/(1-vemin);tau*(1-taumin-tau)/(1-taumin)];

%specify priors
lam_phi=pi/sqrt(4);%variance for gaussian - pi/2 means E(phi)=1
lp_phi=-.5*sum(phi.^2)/lam_phi;    dlp_phi=-(phi.^1)/lam_phi;
a_tau=1;b_tau=2;%beta
lp_tau=(a_tau-1)*log(tau)+(b_tau-1)*log(1-tau);    dlp_tau=(a_tau-1)/tau+(b_tau-1)/(1-tau);
a_ve=2;b_ve=1;%beta
lp_ve=(a_ve-1)*log(ve)+(b_ve-1)*log(1-ve);    dlp_ve=(a_ve-1)/ve+(b_ve-1)/(1-ve);
lp=(sum(lp_phi)+lp_ve+lp_tau);


dlp=[dlp_phi;dlp_ve;dlp_tau];

%construct base covariance matrix
lC0=0;lQ=0;lM=0;lC0j=0;lQj=0;lMj = 0;lQj1=0;lMj1=0;lC0j1=0;
R=1;


z=X(:,1:d)*R;

for i=1:d
    D{i}=abs(z(:,i)*ones(1,T)-ones(T,1)*z(:,i)').^dexp;
    S{i}=(abs(z(:,i)).^dexp)*ones(1,T)+ones(T,1)*abs(z(:,i)').^dexp;
    U{i}=abs(z(:,i)).^dexp;
    lC0=lC0-phi(i)*D{i};
    lQ=lQ-phi(i)*S{i};
    lM=lM-phi(i)*U{i};
    if i== Condj0(1)
        lC0j=lC0j-phi(i)*D{i};
        lQj=lQj-phi(i)*S{i};
        lMj=lMj-phi(i)*U{i};
    end
    if i== Condj1(1)
        lC0j1=lC0j1-phi(i)*D{i};
        S1{i}=(abs(z(:,i)-1).^dexp)*ones(1,T)+ones(T,1)*abs(z(:,i)'-1).^dexp;
        U1{i}=abs(z(:,i)-1).^dexp;
        lQj1=lQj1-phi(i)*S1{i};
        lMj1=lMj1-phi(i)*U1{i};
    end
end

%cond0=0;cond0=1; %
Fcond0=0;

 Cd=tau*(exp(lC0)-cond0*exp(lQ)-condj0*exp(lC0-lC0j+lQj)-condj1*exp(lC0-lC0j1+lQj1));
Md=cond0*exp(lM)*Fcond0+condj0*exp(lMj)*Fcondj0+condj1*exp(lMj1)*Fcondj1;

% %soft cond
% vsoft=0.1;
%  Cd=tau*(exp(lC0)-cond0*exp(lQ)-condj0*exp(lC0-lC0j+lQj)./(1+vsoft./(tau*exp(lC0-lC0j)))-condj1*exp(lC0-lC0j1+lQj1));
% Md=cond0*exp(lM)*Fcond0+condj0*exp(lMj)./(1+vsoft./(tau*exp(lC0-lC0j)))*Fcondj0+condj1*exp(lMj1)*Fcondj1;



mpt=zeros(T,1);

Cdt=Cd;
like=0;
dl=0*dlp;

Id=eye(T);
Sigma=Cd+ve*Id;
dd=det(Sigma);

%chol algorithm from R&W
[L,erp]=chol(Sigma);
a=L\(L'\(Y-Md));
Linv=L\Id;
iKVs=Linv*Linv';
mpt=Md+Cd*a;
Cdt=Cd-Cd*iKVs*Cd;
like=-.5*(Y-Md)'*a-sum(log(diag(L)));

if nargout>1,%calculate gradient
    %a=iKVs*(Y-Md);
    vQ=vec(a*a'-iKVs)';
    for i=1:d
        dC{i}=-D{i}.*exp(lC0)+cond0*S{i}.*exp(lQ);
        dM{i}=-U{i}.*Md;
        dl(i,:)=.5*vQ*vec(dC{i})-dM{i}'*a;
    end
    dC{d+1}=Id;   dl(d+1)=.5*vQ*vec(dC{d+1});
    dC{d+2}=Cd/tau;  dl(d+2)=.5*vQ*vec(dC{d+2});
    
    %J is gradient in parameter space - need gradient in transformed parameters
    J=dl+dlp;
    GradLpost=J.*dpars;
    neglgrad=-GradLpost;
end
if ~isempty(predgrid)
    
    
    %produce mean and variance on grid specified by predgrid
    %uses loop over grid
    [ng,dd]=size(predgrid);
    xs=[predgrid];
    lC0=0;lQ=0;lM=0;lC0j=0;lQj=0;lMj = 0;lQj1=0;lMj1=0;lC0j1=0;
    zg=xs*R;
    for i=1:d
        D{i}=abs(zg(:,i)*ones(1,T)-ones(ng,1)*z(:,i)').^dexp;
        S{i}=(abs(zg(:,i)).^dexp)*ones(1,T)+ones(ng,1)*abs(z(:,i)').^dexp;
        U{i}=abs(zg(:,i)).^dexp;
        lC0=lC0-phi(i)*D{i};
        lQ=lQ-phi(i)*S{i};
        lM=lM-phi(i)*U{i};
        
        if i== Condj0(1)
            lC0j=lC0j-phi(i)*D{i};
            lQj=lQj-phi(i)*S{i};
            lMj=lMj-phi(i)*U{i};
        end
        if i== Condj1(1)
            lC0j1=lC0j1-phi(i)*D{i};
            S1{i}=(abs(zg(:,i)-1).^dexp)*ones(1,T)+ones(ng,1)*abs(z(:,i)'-1).^dexp;
            U1{i}=abs(zg(:,i)-1).^dexp;
            lQj1=lQj1-phi(i)*S1{i};
            lMj1=lMj1-phi(i)*U1{i};
        end
        
    end
    Cs=tau*(exp(lC0)-cond0*exp(lQ)-condj0*exp(lC0-lC0j+lQj)-condj1*exp(lC0-lC0j1+lQj1));
    Ms=cond0*exp(lM)*Fcond0+condj0*exp(lMj)*Fcondj0+condj1*exp(lMj1)*Fcondj1;
    
    
    mps=Ms+Cs*a;
    for i=1:ng,
        Cst(i,:)=tau-Cs(i,:)*iKVs*Cs(i,:)';
    end
    out.pred=mps;
    out.var=Cst+ve*(~noVe==1);
    
    for i=1:d
        Dtilde{i} = (zg(:,i)*ones(1,T)-ones(ng,1)*z(:,i)')';
        grad(:,i)=-2*phi(i)*sum(Dtilde{i}.*(Cs'.*a), 1); % not working for now
    end
    out.grad=grad;
else
    
    %iKVs,a,Cdt,like
    dataGP.iKVs = iKVs;
    dataGP.a = a;
    dataGP.Cdt = Cdt;
    dataGP.like = like;
end

lpost=like+lp;
neglpost=-lpost;

lnL_LOO=.5*sum(log(diag(iKVs)))-.5*sum(a.^2./diag(iKVs));

out.mean=mpt;
out.cov=Cdt;
out.LOO=lnL_LOO;



