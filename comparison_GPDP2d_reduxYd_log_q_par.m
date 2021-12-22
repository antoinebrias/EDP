% %this script compares GP DP and straight DP for two-state models
% %embedding dim of 1 to 2 is assumed - though most models seem to have E>2
% % this code is largely the same as reduxYd but the GP and DP parts are
% carried out in logs instead of the original scale

% %models are f(x,p,v,u)
%1 %juv dd growth, p=[Rmax Gmax Sjuv Sadult sel_juv] 
%2 %two locations, p=[Rmax1 Rmax2 m1 m2 sel_loc2] 
%3 %Ginzburg,      p=[Rmax M a] 
%4 %seas. Ricker,  p=[Rmax amplitude, 2*pi*freq ] 
%5 %Nicholson-bailey,   p=[Rmax lam m 0 bycatch] (N-B with dd growth in host)
%6 %contemporary evolution  p=[Rmax b(exponent) a(sensitivity of het) sel_heterozygote sel_recessive] 
%7 %standard 1-d ricker p=[Rmax]
%8 %delay Hassel p=[Rmax exp 0 0]
%9 %delay logistic p=[Rmax 0 0 0]- no noise for this one.
%10 2-species shepherd p=[Rmax2 a12 a21 b sel] %rmax1 = 1.5 *rmax2

%this version does the gp with (Yt,Yt-1) as inputs (not S as in_reduxB,C)

%version_q allows for bycatch





close all
clear
rng('default');
rng(1988);



% first_iter = 1;
idx_loop =0;
iter_max = 100;
parset_length = 3;
harvestfun_length = 1;
bycatch_length = 1;
model_length = 4;
vset_length = 3;
Tset_length = 3;

resTmp = zeros(iter_max*parset_length*harvestfun_length*bycatch_length*model_length*vset_length*Tset_length,4);

for parset=1:3
for harvestfun=3    
    
%% set up models
switch parset
    case 1% %chaotic/quasicycle regime
    p=[35 .9 .5 .1 0; %invariant loop - quasi cycle; need to check ref (
       15 18 .65 .75 0;%chaos 
       3.5 10 5 0 0; %chaos - looks like Ricker
       log(7) .4 2*pi/4 0 0;%chaos
       7 2 .9 0 0;%invariant loop  
       5 4 .9 .75 0;%chaos, 'strange' attractor
       15 0 0 0 0;
       4 3 0 0 0;%invariant loop
       2.1 0 0 0 0;
       3 .2 .1 4 0];

    case 2 % %limit cycle
    p=[34 .9 .5 .1 0.2;%stable 4-cycle
       13 18 .65 .75 .1;%stable 4-cycle
       4 10 5 0 0;% stable 4-cycle
       log(7) .4 2*pi/3 0 0;%stable 4-cycle
       8 2 .9 0 .1;%stable 5-cycle
       5 3.35 .9  .65 .3;%stable 4-cycle
       12 0 0 0 0;%stable 4-cycle
       4 2 0 0 0;%stable 6-cycle
       2.19 0 0 0 0;
       2.4 .2 .1 4 0];

    case 3, %fixed point
    p=[24 .9 .5 .1 0.2;%fixed point
       7 18 .65 .75 .1;%fixed point
       1.75 10 5 0 0;%fixed point
       log(7) .4 2*pi/8 0 0;%two-cycle.  fixed point would be stupid for seasonal model
       5 2 .9 0 .1; %8.2 .9 -1 .9 .1
       5 2.5 .9  .65 .3;
       15 0 0 0 0;
       4 1.75 0 0 0;
       1.54 0 0 0 0;
       1.6 .2 .1 4 0];
end
x0={[1 10];
    [log(16) log(15)];
    [5 (1-1/2)^(1/2)];
    [log(10) 0.1];
    [1 1]; 
    [2 4 2];
    [2];
    [1 .8];
    [.2 .75];
    [.55 .88]};


%define model
vlist=[0 .1 .2];
Tlist=[30 50 100];
for bycatch=1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!:6, Imprecise selectivty two location scenario and competition scenario only.
    
for model=[2 3 4 10]
    
%      if first_iter 
%         model = 4;
%     end
    
    
for vset=1:3
    
%      if first_iter 
%         vset = 3;
%     end
    
     v=vlist(vset);
for Tset=1:3
    
%      if first_iter 
%         Tset = 2;
%     end
    
     T=Tlist(Tset);


pm=p(model,:);
pm(5)=(bycatch-1)/10;
%v=0.0;
%T=200;
Tf=100;
dua=.01;dub=0.005;
dusd=0.00;
disc=.95;
nx=100;nu=50;
duphi=.98; duv=5;
initsd=.05;

 U=0.0001*ones(1,T);
switch harvestfun
    case 0, %constant
        uinc=@(Yield,U,t) 0;
    case 1, %Fish-o-Stat
        uinc=@(Yield,U,t) (dua*(Yield(t)>=Yield(t-1))*(1-U(t))-dub*(Yield(t)<Yield(t-1))*U(t))*3*(t/T)^2;
    case 2, %stepped
        uinc=@(Yield,U,t) 5*.5/T*(mod(t,5)==0);
    case 3, %AR to u=0.5
        uinc=@(Yield,U,t) 1/(1+((1-U(t))/U(t))^duphi*exp(-(1-duphi)*randn(1)*duv))-U(t);
end
% u=.2*rand(1,T);
% u(1)=0.001;u(2)=.01;
% du=0;

idx_loop = idx_loop+1;
res_par_Tmp = zeros(iter_max,4);

% for iter=1:iter_max   
 parfor iter=1:iter_max   
%       if first_iter 
%         Tset = 94;
%     end
[M,H]=nph_models(model,pm,v);
    
% pm(5)=1

[parset harvestfun model v T iter];

%% steady state yield
Ts=500;
us=linspace(0,.99,100);
xs=zeros(2+(model==6)-(model==7),Ts);
xs(:,1)=x0{model}';
harvs=[];yavg=[];
for j=1:length(us);
    for t=1:Ts
        xs(:,t+1)=M(xs(:,t),us(j));
        harvs(t)=H(xs(:,t),us(j));
    end
    %figure(7);plot([0:Ts],xs(1,:));pause(.1)
    yavg(j)=mean(harvs);
end
ys_opt=max(yavg);us_opt=us(yavg==ys_opt);
 % figure(1);subplot(2,2,1);plot(us,yavg);
    for t=1:Ts
        xs(:,t+1)=M(xs(:,t),us_opt);
    end
% subplot(2,1,2);plot([0:Ts],xs(1,:))
% 
%% simulate data
%this section generates a simulated time series and evaluates a GP forecast
newU=@(u,du) max(0.001,min(1-.001,u+du));%constant F


xt=zeros(2+(model==6)-(model==7),T);
xt(:,1)=x0{model}'*(1+initsd*2*(rand(1)-.5));
Yield=zeros(1,T);
h=zeros(1,T);u = U;

%simulate time series
for t=1:T
    xt(:,t+1)=M(xt(:,t),u(t));
    Yield(t)=H(xt(:,t),u(t));
    h(t)=Yield(t)*(1/u(t)-1);%index of escapement, like X(t)-Y(t) in scalar model
    if t>2
%         du=dua*(Yield(t)-Yield(t-1))/(2*dua+Yield(t)+Yield(t-1))*(1-u(t));
%         du=(dua*(Yield(t)>Yield(t-1))*(1-u(t))-dub*(Yield(t)<Yield(t-1))*u(t))*3*(t/T)^2;
        du=uinc(Yield,u,t)*(1+dusd*2*(rand(1)-.5));
        u(t+1)=newU(u(t),du);
    end
end
% figure(2); subplot(2,2,1);plot([0:T],xt(:,1:T+1),'b',[0:T-1],Yield(1,1:T),'r');
% subplot(2,2,2);plot([0:T],u(1:T+1))
% subplot(2,2,3);plot3(xt(1,1:T-1)-Yield(1,1:T-1),xt(1,2:T)-Yield(1,2:T),xt(1,3:T+1),'.')
% subplot(2,2,4);plot3(h(1,1:T-2),h(1,2:T-1),xt(1,3:T),'.')


% 
% figure(3);subplot(2,2,1);plot([0:T],xt);subplot(2,2,2);hist(xt(1,T-200:T),100);
% subplot(2,2,3);plot([0:T-1],Yield(1,1:T),'r');subplot(2,2,4);plot([0:T-1],u(1,1:T),'r')
% 
% figure(4);plot(xt(1,T-100:T),xt(1,T-100+1:T+1),'.');
% figure(2);plot3(xt(1,1:T-1),xt(1,2:T),xt(1,3:T+1),'.');



% pause
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up grids
%lb=min(0.0001,min(h));ub=max(h)+2*std(h);
lb=log(.75*min(h));ub=log(max(h)*1.25);
hg=linspace(lb,ub,nx)';
dh=hg(2)-hg(1);
lowu=.0001;hiu=1-.0001;
ug=linspace(lowu,hiu,nu);ugt=ug';
du=ug(2)-ug(1);
lnhgu=log(exp(hg)*(1./(1-ug)))';
%lnhgu(nu,:)=-20*(hg'>lb);
%hgtile=ones(nu,1)*hg';


%% parametric model fit
Btp1=(h(2:T)+Yield(2:T))';Bt=h(1:T-1)';
%rx=[ones(T-1,1) Bt log(Bt)];
%ry=log(Btp1);
rx=[ones(T-1,1) Bt];
ry=log(Btp1./Bt);
pt_pars=inv(rx'*rx)*rx'*ry;
%Btp1_pred=exp(rx*pt_pars);
Btp1_pred=exp(log(Bt)+rx*pt_pars);
pt_err=(ry-rx*pt_pars)'*(ry-rx*pt_pars)/(T-1);
covpars=inv(rx'*rx);
%   figure(28);plot(Bt,Btp1,'.',Bt, Btp1_pred,'r.')
% figure(29);plot3(Bt(1:end-1),Bt(2:end),Btp1(2:end),'.')


%% Bayesian Parametric transition matrix
msurf_par=[];vsurf_par=[];
P_par=[];P0_par=[];
for i=1:nx,%lognormal transitions
    rxi=[1 exp(hg(i)) (hg(i))];
    pt_pars(3)=1;covpars(3,3)=0;
    msurf_par(i)=rxi*pt_pars;
%    vsurf_par(i)=pt_err*([1 hg(i) log(hg(i))]*covpars*[1 hg(i) log(hg(i))]'+1);
    vsurf_par(i)=max(dh/2,pt_err*(rxi*covpars*rxi'+1));%this overcomes the problem of 0 probability when the error is smaller than the grid spacing
    ff=exp(-.5*(lnhgu-msurf_par(i)).^2/vsurf_par(i));
    ff(isnan(ff))=0;
    
    z1=(lnhgu(:,1)-msurf_par(i))/vsurf_par(i)^.5;
    ff(:,1)=ff(:,1)+normcdf(z1)/dh;
    
    z2=(lnhgu(:,nx)-msurf_par(i))/vsurf_par(i)^.5;
    ff(:,nx)=ff(:,nx)+(1-normcdf(z2))/dh;
    
    P_par{i}=diag(1./sum(ff,2))*ff;
    P0_par(i,:)=P_par{i}(1,:);
end;
Fbar_par=exp(msurf_par+0.5*vsurf_par);

% figure(31);plot(log(Bt),Btp1,'k.',hg,Fbar_par,'k')
% figure(32);plot(hg,msurf_par);  
% surf(hg,hg,P0_par)


%% parametric SDP
vs=exp(hg)'-disc*Fbar_par;
[c,ind]=min(vs);zs=hg(ind);
V_par=hg-c;
uopt_par=max(0,1-exp(zs-hg));
%V_par=Fbar_par'/(1-disc);
%V_par(1)=0;
count=0;deltau=1;
%iterate 
unew_par=[];
while (deltau>=du)&(count<1000)  
%figure(999);subplot(2,1,1);plot(hg,V_par);subplot(2,1,2);plot(hg,uopt_par);%pause
    for i=1:nx
        [Vstar,ind]=max(Fbar_par(i)*ugt+disc*P_par{i}*V_par);
        V_par(i)=Vstar;
        unew_par(i)=ug(ind);
    end
    deltau=max(abs(uopt_par(:)-unew_par(:)));
    uopt_par=unew_par;
    count=count+1;
end
%  figure(34);plot(hg,uopt_par); 
% figure(35);plot(hg,V_par); 
%  
%% Evaluate parametric policy 
htnew_par=zeros(1,T+Tf);Yield_par=zeros(1,T+Tf);u_par=zeros(1,T+Tf);xt_par=zeros(2+(model==6)-(model==7),T+Tf);
htnew_par(1:T)=log(h);
u_par(1:T)=u(1:T);Yield_par(1:T)=Yield;xt_par(:,1:T+1)=xt;
for t=T+1:T+Tf
    inext=max(1,min(nx,floor((htnew_par(t-1)-lb)/dh)+1));
    u_par(t)=uopt_par(inext);
    Yield_par(t)=H(xt_par(:,t),u_par(t));
    htnew_par(t)=log((Yield_par(t)*(1/u_par(t)-1)));
    xt_par(:,t+1)=M(xt_par(:,t),u_par(t));
end


% figure(37);subplot(3,1,1);plot([0:T+Tf],xt_par(1,:),'b');
%            subplot(3,1,2);plot([1:T+Tf],u_par,'b',[0 T+Tf],us_opt*[1 1],'r');
%            subplot(3,1,3);plot([1:T+Tf],Yield_par,'b',[0 T+Tf],ys_opt*[1 1],'r',[T+1:T+Tf],mean(Yield_par(T+1:T+Tf))*ones(Tf,1),'g');
    
% pause
%%           
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% evaluating the embedding dimension;
% %yield and effort as seperate inputs?
% X=([Yield(2:T-1)' u(2:T-1)' Yield(1:T-2)' u(1:T-2)']);
% Y=log(h(3:T)+Yield(3:T))';
% for d=1:4
%     lpost=@(p) GP4DP_Embed2b_notheta2(p,X(2:end,1:d),Y(2:end),1,[],d);
%     LenScale=.1*ones(1,d);
%     ve=.000000001; tau=10;
%     pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';%theta fixed
% %     figure(17)
%     [pf{d},nlogl(d)]=fmingrad_Rprop(lpost,pars);
% end
% nlogl

% % embedding dim from 1 to 4
% X=([h(4:T-1)' h(3:T-2)' h(2:T-3)' h(1:T-4)']);%3 or 4 lags of h?
% Y=log(h(5:T)+Yield(5:T))';
% mY=min(Y);sdY=max(Y)-mY;
% Y=(Y-mY)/sdY;
% for d=1:4
%     lpost=@(p) GP4DP_Embed2b_fit(p,X(1:end,1:d),Y(1:end),0,[],d,0);
%     LenScale=.1*ones(1,d);
%     ve=.000000001; tau=10;
%     pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';%theta fixed
% %     figure(17)
%     [pf{d},nlogl(d)]=fmingrad_Rprop(lpost,pars);
% end
% modelE(model,:)=nlogl;

    

%for fitting model (just E=1 or 2)
X=log([h(2:T-1)' h(1:T-2)']);%3 or 4 lags of h?


Y=log(h(3:T)+Yield(3:T))';

% [X,Y]

mY=mean(Y);sdY=std(Y);
Y=(Y-mY)/sdY;
% [X,Y]
% size(X)
pf=[];nlogl=[];
for d=1:2
    lpost=@(p) GP4DP_Embed2b_fit(p,X(1:end,1:d),Y(1:end),0,[],d,0,[0 0]);
%     [X(1:end,1:d),Y(1:end)]
    LenScale=.05*ones(1,d);
    ve=.000000001; tau=10;
    pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';%theta fixed
%     figure(17)
    [pf{d},nlogl(d)]=fmingrad_Rprop(lpost,pars);
end
bestE=2;
if nlogl(2)>(nlogl(1)-2), bestE=1;end 
% 
if nlogl(2)>(nlogl(1)-2), 
        LenScale=.05*ones(1,1);
    ve=.000000001; tau=10;
    pars1=[log(LenScale) log(ve/(1-ve)) log(tau)]';%theta fixed
    [dummy,dummy2,out]=GP4DP_Embed2b_fit(pars1,X(:,1),Y,0,[X(:,1)],1,0,[]);
    mp1=out.pred;
         
    LenScale=.05*ones(1,2);
    ve=.000000001; tau=10;
    pars2=[log(LenScale) log(ve/(1-ve)) log(tau)]';%theta fixed
    [dummy,dummy2,out]=GP4DP_Embed2b_fit(pars2,X(:,1:2),Y,0,[X(:,1:2)],2,0,[]);
    mp2=out.pred;

%     figure(999); 
%     subplot(2,2,1);plot([3:T],exp(Y),'.',[3:T],exp(mp1),'ko',[3:T],exp(mp2),'ro')
%     subplot(2,2,2);plot3(X(:,1),X(:,2),exp(Y),'.',X(:,1),X(:,2),exp(mp1),'ko')
%     subplot(2,2,4);plot3(X(:,1),X(:,2),exp(Y),'.',X(:,1),X(:,2),exp(mp2),'ko')
% %     keyboard
end


%% determine GP
%initialize parameters
EmbedDim=bestE;bestRot=[];
if bestE==1, bestRot=0;end
LenScale=.05*ones(1,EmbedDim);
ve=.000000001;
tau=10;
pars=[log(LenScale) log(ve/(1-ve)) log(tau)]';

%create data matrix and rescale
X=log([h(2:T-1)' h(1:T-2)']);
Y=log(h(3:T)+Yield(3:T))';
mY=mean(Y);%min(Y);
sdY=std(Y);%max(Y)-min(Y);
Y=(Y-mY)/sdY;

%n-d embedding on original axes   
lpost=@(p) GP4DP_Embed2b_fit(p,X(:,1:EmbedDim),Y,0,[],EmbedDim,0,[]);
[pfit,nll]=fmingrad_Rprop(lpost,pars);
[nll,ngl,out]=GP4DP_Embed2b_fit(pfit,X(:,1:EmbedDim),Y,0,[X(:,1:EmbedDim)],EmbedDim,0,[]);
Cd=out.var; 
mp=out.pred;
    
% figure(3); 
% subplot(2,1,1);plot([3:T],exp(Y),'.',[3:T],exp(mp),'ko')
% subplot(2,1,2);plot3(X(:,1),X(:,2),exp(Y),'.',X(:,1),X(:,2),exp(mp),'ko')
% figure(4);plot(X(:,1),exp(Y),'.',X(:,1),exp(mp),'ko')


% %2-d embedding on 45 degree rotation
if bestE==2,
    lpost=@(p) GP4DP_Embed2b_fit(p,X,Y,0,[],EmbedDim,(pi/4),[]);
    [pfit2,nll2]=fmingrad_Rprop(lpost,pars);
    [nll2,ngl,out]=GP4DP_Embed2b_fit(pfit2,X,Y,0,[X],EmbedDim,(pi/4),[]);
    Cd=out.var; 
    mp=out.pred;
    bestRot=pi/4*(nll2<nll);pfit=pfit+(pfit2-pfit)*(nll2<nll);

%     figure(5); 
%     subplot(2,1,1);plot([3:T],exp(Y),'.',[3:T],exp(mp),'ko')
%     subplot(2,1,2);plot3(X(:,1),X(:,2),exp(Y),'.',X(:,1),X(:,2),exp(mp),'ko')
%     figure(6);plot(X(:,1),exp(Y),'.',X(:,1),exp(mp),'ko')
end

% %n-d embedding on best axes: effects of conditioning f(0,0)   
best_cond_pars=[0 0];
% condfactor=[0 .25 .5 .75 1]';
% condpars=[condfactor -2./(condfactor+(condfactor==0))];
% for i=1:length(condfactor)                    
%     lpost=@(p) GP4DP_Embed2b_fit(p,X(:,1:EmbedDim),Y,0,[],EmbedDim,bestRot,condpars(i,:));
%     [pfit_cond,nll_cond]=fmingrad_Rprop(lpost,pars);
%     [nll_cond,ngl,out_cond]=GP4DP_Embed2b_fit(pfit_cond,X,Y,0,[X],EmbedDim,bestRot,condpars(i,:));
%     Cond_LOO(i)=out_cond.LOO;
% end
% ind=find(Cond_LOO==max(Cond_LOO))
% best_cond_pars=condpars(ind,:);


% figure(333); 
% subplot(2,1,1);plot([3:T],exp(Y),'.',[3:T],exp(mp),'ko')
% subplot(2,1,2);plot3(X(:,1),X(:,2),exp(Y),'.',X(:,1),X(:,2),exp(mp),'ko')
% figure(4);plot(X(:,1),exp(Y),'.',X(:,1),exp(mp),'ko')


% % % 
%% Projection matrix from GP4DP_Embed2b_fit

%create transition matrix from GP, for each u
%get predictions by matrix
[hg1,hg2]=meshgrid(hg,hg);
vg1=vec(hg1');vg2=vec(hg2');

if EmbedDim==1,
    [nll,ngl,out]=GP4DP_Embed2b_fit(pfit,X(:,1:EmbedDim),Y,0,hg,EmbedDim,bestRot,best_cond_pars);
    mh=out.pred*sdY+mY;
    vh=out.var*sdY^2;
    vsurf=vh*ones(1,nx);
    msurf=mh*ones(1,nx);
else
    [nll,ngl,out]=GP4DP_Embed2b_fit(pfit,X(:,1:EmbedDim),Y,0,[vg1 vg2],EmbedDim,bestRot,best_cond_pars);
    mh=out.pred*sdY+mY;
    vh=out.var*sdY^2;
    vsurf=reshape(vh,nx,nx);%make sure that v includes ve
    msurf=reshape(mh,nx,nx);
end



% msurf=reshape(log(pm(1)*vg1.*(1-vg2)),nx,nx);
% vsurf=.00001*ones(nx);

P=[];Fbar_num=[];
for i=1:nx,for j=1:nx,%normal transitions
%    ff=exp(-.5*((lnhgu-mY)/sdY-msurf(i,j)).^2/vsurf(i,j))./hgtile;
    ff=exp(-.5*(lnhgu-msurf(i,j)).^2/vsurf(i,j));
    ff(isnan(ff))=0;
    
    z1=(lnhgu(:,1)-msurf(i,j))/vsurf(i,j)^.5;
    ff(:,1)=ff(:,1)+normcdf(z1)/dh;
    
    z2=(lnhgu(:,nx)-msurf(i,j))/vsurf(i,j)^.5;
    ff(:,nx)=ff(:,nx)+(1-normcdf(z2))/dh;
   
    sf=sum(ff,2);
    P{i,j}=diag(1./(sf+(sf==0)))*ff;
% figure(999);surf(ug,hg,P{i,j}');pause(.01)
%shading flat;hold on;plot(ug,exp(msurf(i,j))*(1-ug),'w');hold off;   
    Fbar_num(i,j)=P{i,j}(1,:)*exp(hg);
end;end
Fbar=Fbar_num;
% figure(11);surf(hg,hg,vsurf);hold on;plot(X(:,1),X(:,2),'ko');hold off 
% pfit
% figure(12);surf(hg,hg,msurf);hold on;plot3(X(:,2),X(:,1),Y*sdY+mY,'ko');hold off    
% figure(120);surf(exp(hg),exp(hg),exp(msurf)');hold on;plot3(exp(X(:,1)),exp(X(:,2)),exp(Y*sdY+mY)','ro');hold off
% exp(pfit(1:2))



% % 
%% determine value function
deltau=1; count=0;
if bestE==1,
    V=V_par;Vnew=V;
    uopt=zeros(nx,1);
    unew=uopt;

%     figure(14);
%     subplot(2,2,1);plot(hg,Fbar);hold on;plot(h(2:T),h(1:T-1),'w.');hold off;
%     subplot(2,2,4);plot(hg,V);
    %iterate 
    while (deltau>=du)&(count<5000)  
        for i=1:nx, 
           [Vstar,ind]=max(Fbar(i)*ugt+disc*P{i,1}*V);
            Vnew(i)=Vstar;
            unew(i)=ug(ind);
        end
        deltau=max(abs(uopt(:)-unew(:)));
        uopt=unew;
        V=Vnew;V(1,1)=0;
%         figure(14);subplot(2,2,3);plot(hg,uopt);subplot(2,2,4);plot(hg,V); 
        count=count+1;
    end
    uopt=uopt*ones(1,nx);
else, %DP for 2-d embedding
    %initialize value function
    V=V_par*ones(1,nx);Vnew=V;
    uopt=zeros(nx,nx);
    unew=uopt;

%     figure(14);
%     subplot(2,2,1);pcolor(hg,hg,Fbar);shading flat;colorbar;hold on;plot(log(h(2:T)),log(h(1:T-1)),'w.');hold off;
%     subplot(2,2,4);pcolor(hg,hg,V);shading flat;colorbar;pause(.1)

    %iterate 
    while (deltau>=du)&(count<1000)  
        for i=1:nx, 
            Vi=V(:,i);
            for j=1:nx
                [Vstar,ind]=max(Fbar(i,j)*ugt+disc*P{i,j}*Vi);
                Vnew(i,j)=Vstar;
                unew(i,j)=ug(ind);
            end
        end
        deltau=max(abs(uopt(:)-unew(:)));
        uopt=unew;
        V=Vnew;
%          figure(14);subplot(2,2,3);pcolor(hg,hg,uopt);colorbar;shading flat; ;subplot(2,2,4);pcolor(hg,hg,V);colorbar;shading flat;pause(.1) 
        count=count+1;
        
        
%         figure(555);
%              [x1_uopt,x2_uopt]=meshgrid(exp(hg),exp(hg));
%         surf(x1_uopt,x2_uopt,uopt);
%         pause
    end
%       figure(14);subplot(2,2,3);pcolor(hg,hg,uopt);colorbar;shading flat;subplot(2,2,4);pcolor(hg,hg,V);colorbar;shading flat; 
%  
end% 

%% Evaluate DP Policy
ht_GP=zeros(1,T+Tf);Yield_GP=zeros(1,T+Tf);u_GP=zeros(1,T+Tf);xt_GP=zeros(2+(model==6)-(model==7),T+Tf);
ht_GP(1:T)=log(h);
Yield_GP(1:T)=Yield;u_GP(1:T)=u(1:T);xt_GP(:,1:T+1)=xt;
for t=T+1:T+Tf
    inext=max(1,min(nx,floor((ht_GP(t-1)-lb)/dh)+1));
    jnext=max(1,min(nx,floor((ht_GP(t-2)-lb)/dh)+1));
    u_GP(t)=uopt(inext,jnext);
    Yield_GP(t)=H(xt_GP(:,t),u_GP(t));
    ht_GP(t)=log(Yield_GP(t)*(1/u_GP(t)-1));
    xt_GP(:,t+1)=M(xt_GP(:,t),u_GP(t));
end
% figure(550);
%         [x1_uopt,x2_uopt]=meshgrid(exp(hg),exp(hg));
%         surf(x1_uopt,x2_uopt,uopt);
% figure(7); subplot(3,1,1);plot([0:T+Tf],xt_GP(1,:),'b');
%            subplot(3,1,2);plot([1:T+Tf],u_GP,'b',[0 T+Tf],us_opt*[1 1],'r');
%            subplot(3,1,3);plot([1:T+Tf],Yield_GP,'b',[0 T+Tf],ys_opt*[1 1],'r',[T+1:T+Tf],mean(Yield_GP(T+1:T+Tf))*ones(Tf,1),'k',[T+1:T+Tf],mean(Yield_par(T+1:T+Tf))*ones(Tf,1),'g');
%     
%%
% 
%
% 444
%% save stuff

%if model<7, pause;end
% relY(model,:)=[ys_opt mean(Yield_GP(T+1:T+Tf)) mean(Yield_par(T+1:T+Tf)) mean(Yield_GP(T+1:T+Tf))/mean(Yield_par(T+1:T+Tf))];
% Cond_LOO_set(model,:)=Cond_LOO;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TD-EDP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X=log([h(2:T-1)' ]);
Y=log(h(3:T)+Yield(3:T))';

% x_dat = exp(X);%data.xt(1,:)';
% u_dat = (Yield(3:T)./(h(3:T)+Yield(3:T)))';

x_dat = xt(1,:)';
u_dat = u(1,:)';

% default EDM structure
[mdlstruct] = init_mdl_struct(1); % the input is the number of available species data

% data generation
mdlstruct.data = x_dat;
mdlstruct.control_data = u_dat ;

%different kinds of control exist, we need to specify which one we use here
mdlstruct.control_type =  'rate';
mdlstruct.control_param =  [1 ]; % both  dimensions are fully controlled

% % % data display
% figure
% data_show(mdlstruct)


for d=1:2
% we specify the number of lags used here
% here, we specify the number of time lags used for each variables
% 0 means that the variable is not used in the prediction
mdlstruct.n_lags = [d];

% % GP regression
% disp('*****************************************************')
% disp('GP fitting')
% disp('*****************************************************')
% this step can be skipped if a model already exists, otherwise a GP is
% fitted to the data
% fills the field mdlstruct.model with the fitted gp
% fills the field mdlstruct.model_stats with statistics computed during the
% GP regression

% default gp structure
mdlstruct.gp = init_gp();

mdlstruct.gp.is_log = 0;
mdlstruct.gp.cond_0=0;

% %%%%%%%%%%%%%%%!!!!!
 mdlstruct.gp.is_value_function=1;

% we now fit the gp with the data in mdlstruct
mdlstruct = fit_gp(mdlstruct);

nll(d)=mdlstruct.gp.nll;
end
bestEgp=2;
if nll(2)>(nll(1)-2), bestEgp=1;end 

mdlstruct.n_lags = [bestEgp];
mdlstruct.gp = init_gp();
mdlstruct.gp.is_log =0;
mdlstruct.gp.cond_0=0;
mdlstruct = fit_gp(mdlstruct);



mdlstruct.value_function_type = 'gp'; %'gp'or 'nn'
%using logarithmic value function
mdlstruct.is_value_function_log = 0; %not working rn


% % gp display
%  gp_show(mdlstruct)


% % Optimal control problem
% disp('*****************************************************')
% disp('  Optimal control problem structure creation')
% disp('*****************************************************')

optstruct=[];

% Discount factor
optstruct.discount_factor=disc;

% Objective weights
optstruct.weights = [1];


% Reward function
% The instantaneous reward function is a handle function.
% The input are the current variables, the controls and the weight for each
% dimension
% The output are the total weighted reward, the unweighted reward on each
% dimension, and the weighted reward on each dimension

%>>> Example 1 <<<<
% use of an existing function
optstruct.reward = @(x,u,w) reward_harvest_all(x,u,w);

% % Temporal difference learning
% disp('*****************************************************')
% disp('EDP using Temporal difference learning')
% disp('*****************************************************')
% default  td structure
dpstruct=init_td(mdlstruct);


%% using 1sp result to build the initial value function
 dpstruct.value_function_init = fit_value_function([],[(hg)],V_par.*0,mdlstruct);


%  run the td algorithm, using the gp regression
 [dpstruct] = td_learning_v2(mdlstruct,optstruct,dpstruct,'gp');

% %  % control maps
% %  disp('Control maps given by the EDP algorithm')
% td_show_v2(mdlstruct,optstruct,dpstruct)

% controlled trajectories
t_max = 100; % time length
n_traj = 1; % number of trajectories
start_x = 'last'; % starting point
is_disp = 0; % display the trajectories result
% %  disp('Controlled trajectories following the EDP policy')
%  traj = sim_trajectories_v2(mdlstruct,optstruct,dpstruct,t_max,n_traj,start_x,is_disp);



%  [opt_control_tmp,unweighted_opt_value,weighted_opt_value]=td_policy(X{t-1},optstruct,dpstruct,mdlstruct);

 %% Evaluate DP Policy
ht_TD_GP=zeros(1,T+Tf);Yield_TD_GP=zeros(1,T+Tf);u_TD_GP=zeros(1,T+Tf);xt_TD_GP=zeros(2+(model==6)-(model==7),T+Tf);
ht_TD_GP(1:T)=log(h);
Yield_TD_GP(1:T)=Yield;u_TD_GP(1:T)=u(1:T);xt_TD_GP(:,1:T+1)=xt;
for t=T+1:T+Tf
%     inext=max(1,min(nx,floor((ht_TD_GP(t-1)-lb)/dh)+1));
%     jnext=max(1,min(nx,floor((ht_TD_GP(t-2)-lb)/dh)+1));
    
     [opt_control_tmp,unweighted_opt_value,weighted_opt_value]=td_policy_v2(xt_TD_GP(1,t-1:-1:t-2).*(1-u_TD_GP(1,t-1:-1:t-2)),optstruct,dpstruct,mdlstruct);
    
    u_TD_GP(t)=opt_control_tmp;%uopt(inext,jnext);
    Yield_TD_GP(t)=H(xt_TD_GP(:,t),u_TD_GP(t));
    ht_TD_GP(t)=log(Yield_TD_GP(t)*(1/u_TD_GP(t)-1));
    xt_TD_GP(:,t+1)=M(xt_TD_GP(:,t),u_TD_GP(t));
end
% figure(77); subplot(3,1,1);plot([0:T+Tf],xt_TD_GP(1,:),'b',[0:T+Tf],xt_GP(1,:),'c');
%            subplot(3,1,2);plot([1:T+Tf],u_TD_GP,'b',[1:T+Tf],u_GP,'c',[0 T+Tf],us_opt*[1 1],'r');
%            subplot(3,1,3);plot([1:T+Tf],Yield_TD_GP,'b',[0 T+Tf],ys_opt*[1 1],'r',[T+1:T+Tf],mean(Yield_TD_GP(T+1:T+Tf))*ones(Tf,1),'k',[T+1:T+Tf],mean(Yield_par(T+1:T+Tf))*ones(Tf,1),'g',[T+1:T+Tf],mean(Yield_GP(T+1:T+Tf))*ones(Tf,1),'c');
    
%  [mean(Yield_par(T+1:T+Tf)) mean(Yield_GP(T+1:T+Tf)) mean(Yield_TD_GP(T+1:T+Tf)) ];
 
%  resYield(iter,Tset,vset,model,bycatch,harvestfun,parset,1:3)=[mean(Yield_par(T+1:T+Tf)) mean(Yield_GP(T+1:T+Tf)) mean(Yield_TD_GP(T+1:T+Tf)) ];
 
%  resTmp = [resTmp; [ys_opt mean(Yield_par(T+1:T+Tf)) mean(Yield_GP(T+1:T+Tf)) mean(Yield_TD_GP(T+1:T+Tf)) ] ];
 res_par_Tmp(iter,:) = [ys_opt mean(Yield_par(T+1:T+Tf)) mean(Yield_GP(T+1:T+Tf)) mean(Yield_TD_GP(T+1:T+Tf)) ]
% filename=['C:\Users\Renaud\Documents\MATLAB\Steve\resComparison_test\output2_' num2str(parset) '_' num2str(harvestfun) '_' num2str(bycatch) '_' num2str(100*v) '_' num2str(T) '_' num2str(iter)];
% save(filename,'ys_opt','us_opt','xt_GP','xt_par','xt_TD_GP','u_TD_GP','Yield_TD_GP','Yield_GP','Yield_par','pfit','u_GP','u_par','uopt','uopt_par','model','pm','pt_pars','pt_err','bestE')

%    if first_iter 
%         first_iter = 0;
%     end



 end
  resTmp((idx_loop-1)*iter_max+1:iter_max,:) = res_par_Tmp;
 
 end;end;
model
end

end;
harvestfun
end
parset
end

save('december_res_par.mat')

% modelE
% apxE=1+sum(diff(-modelE,[],2)>4,2);
% [apxE relY(:,4)]

 
% %            