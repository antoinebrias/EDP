function [mu,v] = dyn_ricker(x,u,is_det,param)
% DYN_RICKER  Multi-species Ricker model.
%   [MU,V] = DYN_RICKER(X,U,IS_DET,PARAM) computes the next state of a 
%   system in a current state X following a Ricker dynamics. U is the
%   fraction of X catched. IS_DET is a boolean indicating if the result is
%   the deterministic mean or a sample. PARAM is a cell array containing the
%   parameters of the Ricker model. PARAM{1} is
%   the standard deviation of the log-normal noise applied to the model.
%   PARAM{2} is the interaction matrix and PARAM{3} is the growth rate array.

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
