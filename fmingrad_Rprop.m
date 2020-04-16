function [xopt,fopt,gradopt]=fmingrad_Rprop(fun,xinit)
% FMINGRAD_RPROP  R-prop algorithm optimizing the GP hyperparameters.
%   [XOPT,FOPT,GRADOPT]=FMINGRAD_RPROP(FUN,XINIT)
%     FUN is a handle to a function that has grad as optional second output
%     this uses the sign of the gradient to determine descent directions 
%     and an adaptive step size - supposedly smoother convergence than
%     conjugate gradients for GP optimization. XINIT is the initial point.
    p=length(xinit);
    
    %optimization parameters for Rprop
    Delta0 = 0.1*ones(p,1);
    Deltamin = 1e-6;
    Deltamax = 50;
    eta_minus = 0.5;eta_minus=eta_minus-1;
    eta_plus = 1.2;eta_plus=eta_plus-1;    
    maxcount=200;
    
    
    %initialize 
    x=xinit;
    [f,g]=fun(xinit);
    s=sqrt(g'*g);
 
    %loop starts here
    count=0;
    del=Delta0;
    df=10;
    while (s>.0001)&(count<100)&(df>.0000001)
        %step 1-move
        xnew=x-sign(g).*del;
        [fnew,gnew]=fun(xnew);
        s=sqrt(gnew'*gnew);
        df=abs(fnew/f-1);

        %step 2 - update step size
        gc=g.*gnew;
        del=min(Deltamax,max(Deltamin,del.*(1+eta_plus*(gc>0)+eta_minus*(gc<0))));
        
        x=xnew;
        g=gnew;
        f=fnew;
        count=count+1;
    end

xopt=x;[fopt,gradopt]=fun(x);

    
    