function[neglpost,neglgrad,out]=GP4DP_Embed2b_fit(pars,xd,yd,Fig,predgrid,EmbedDim, theta,condpars);
%This version is for use in 2-d dynamic programming
%assumes xd is scaled [0,1] and yd is centered
%embedding dimension is 1 to 2, covariance function is
%exponential (dexp=1), squared-exponential (dexp=2), or conditional on f(0,0)=0 (cond0=1);
d=EmbedDim;h=1;
T=length(xd);
Y=yd;
dexp=2;
if ~isempty(condpars),%Fcond*cond is E(f(0)), cond is between 0 and 1 given by (1-var(f(0))/tau) 
    cond0=condpars(1);Fcond0=condpars(2); 
else
    cond0=0;Fcond0=0;
end
npars=d+2;

%pars

%pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';

%transform parameters from real line to constrained space
vemin=0.001;taumin=.1;
phi=exp(pars(1:d));%phi=0.1*ones(d,1);
ve=(1-vemin)*exp(pars(d+1))/(1+exp(pars(d+1)))+vemin;
tau=(1-taumin)*exp(pars(d+2))/(1+exp(pars(d+2)))+taumin;
dpars=[phi;ve*(1-vemin-ve)/(1-vemin);tau*(1-taumin-tau)/(1-taumin)];

%specify priors 
    lam_phi=pi/2;%variance for gaussian - pi/2 means E(phi)=1
    lp_phi=-.5*sum(phi.^2)/lam_phi;    dlp_phi=-(phi.^1)/lam_phi;
    a_tau=1;b_tau=2;%beta
    lp_tau=(a_tau-1)*log(tau)+(b_tau-1)*log(1-tau);    dlp_tau=(a_tau-1)/tau+(b_tau-1)/(1-tau);
    a_ve=2;b_ve=1;%beta
    lp_ve=(a_ve-1)*log(ve)+(b_ve-1)*log(1-ve);    dlp_ve=(a_ve-1)/ve+(b_ve-1)/(1-ve);
    lp=(sum(lp_phi)+lp_ve+lp_tau);
    dlp=[dlp_phi;dlp_ve;dlp_tau];
    
%construct base covariance matrix
lC0=0;lQ=0;lM=0;
R=1;
if EmbedDim==2,
    R=[cos(theta) sin(theta);-sin(theta) cos(theta)];
end

z=xd(:,1:d)*R;

for i=1:d
    D{i}=abs(z(:,i)*ones(1,T)-ones(T,1)*z(:,i)').^dexp;
    S{i}=(abs(z(:,i)).^dexp)*ones(1,T)+ones(T,1)*abs(z(:,i)').^dexp;
    U{i}=abs(z(:,i)).^dexp;
    lC0=lC0-phi(i)*D{i};
    lQ=lQ-phi(i)*S{i};
    lM=lM-phi(i)*U{i};
end

Cd=tau*(exp(lC0)-cond0*exp(lQ));
Md=cond0*exp(lM)*Fcond0;
mpt=zeros(T,1);

Cdt=Cd;
like=0;
dl=0*dlp;

Id=eye(T);
Sigma=Cd+ve*Id;
dd=det(Sigma);
% if dd>1e-6, 
%     iKVs=inv(Sigma);
%     logdd=log(dd);
% else
%     [UU,SS,VV]=svd(Sigma);
%     digS=diag(SS);keptS=(digS>1e-6);
%     iSS=diag((digS>1e-6)./digS);
%     iKVs=VV*iSS*UU';
%     logdd=.5*sum(log(digS+(1-keptS)).*keptS);
% end

% %numerically stable svd for inversion
% [UU,SS,VV]=svd(Sigma);
% digS=diag(SS);keptS=(digS>1e-6);
% iSS=diag((digS>1e-6)./digS);
% iKVs=VV*iSS*UU';
% logdd=.5*sum(log(digS+(1-keptS)).*keptS);

% like=-.5*(Y-Md)'*iKVs*(Y-Md)-.5*logdd;
% %if isinf(like),keyboard;end
% mpt=Cd*iKVs*(Y-Md);
% Cdt=Cd-Cd*iKVs*Cd;



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
    lC0=0;lQ=0;lM=0;
    zg=xs*R;
    for i=1:d
        D{i}=abs(zg(:,i)*ones(1,T)-ones(ng,1)*z(:,i)').^dexp;
        S{i}=(abs(zg(:,i)).^dexp)*ones(1,T)+ones(ng,1)*abs(z(:,i)').^dexp;
        U{i}=abs(zg(:,i)).^dexp;
        lC0=lC0-phi(i)*D{i};
        lQ=lQ-phi(i)*S{i};
        lM=lM-phi(i)*U{i};
    end
    Cs=tau*(exp(lC0)-cond0*exp(lQ));
    Ms=cond0*exp(lM)*Fcond0;
%    mps=Ms+Cs*iKVs*(Y-Md);
    mps=Ms+Cs*a;
    
    for i=1:ng,
      Cst(i,:)=tau-Cs(i,:)*iKVs*Cs(i,:)';
    end
    out.pred=mps;
    out.var=Cst+ve;
end
        
lpost=like+lp;
neglpost=-lpost;

lnL_LOO=.5*sum(log(diag(iKVs)))-.5*sum(a.^2./diag(iKVs));

out.mean=mpt;
out.cov=Cdt;
out.LOO=lnL_LOO;

if Fig==1,
    figure;
    subplot(3,1,1);plot([1:T],mpt,'k',[1:T],Y,'b.');
    subplot(3,1,2);plot(mpt,Y,'.');
    if d==1,
        subplot(3,1,3);plot(xd(1:T-1),mpt(1:T-1),'k.',xd(1:T-1),Y(1:T-1),'b.');
    elseif d==2,
        subplot(3,1,3);plot3(xd(1:T-1),xd(2:T),mpt(2:T),'k.',xd(1:T-1),xd(2:T),Y(2:T),'b.');
    end
    pause(.1)
end
