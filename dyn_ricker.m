%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Multi-species Ricker model
%
% Input are: the population x, the control u, the ricker model parameters
% param{1} contains the standard deviation
% param{2} contains the interaction matrix 
% param{3} contains the growth rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mu,v] = dyn_ricker(x,u,is_det,param)
if ~isempty(u)
    s=x.*(1-u);
else
    s = x;
end
n=size(x,2);
m=size(x,1);
mu = zeros(m,n);

sigma= param{1};
A= param{2};
r= param{3};


for i=1:m
    mu(i,:) =  mdlnd(s(i,:),sigma,A,r,is_det);
end

v = repmat(sigma.^2,size(x,1),1);
end


function y = mdlnd(x,sigma,A,r,noNoise)

% y=x.*(b-A*x')';
y=x.*exp(r-(A*x')');

y(y<0) = 0;
% Add uncertainties
y = y.*exp((1-noNoise)*normrnd(0,sigma,size(x,1),size(x,2)));

end
