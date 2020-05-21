function [mu,v] = dyn_migration(x,u,is_det,param)
% DYN_COMPETITION  Migration model.
%   [MU,V] = DYN_COMPETITION(X,U,IS_DET,PARAM) computes the next state of a 
%   system in a current state X following a migration dynamics. U is the
%   fraction of X catched. IS_DET is a boolean indicating if the result is
%   the deterministic mean or a sample. PARAM is a cell array containing the
%   parameters of the migration model. PARAM{1} is
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
m= param{2};
r= param{3};
b= param{4};


mu = s*diag(r).*exp((1-is_det)*normrnd(0,sigma,size(s,1),size(s,2))).*(1+(s*[1 m(1); m(2) 1]).^b).^-1;
v = repmat(sigma.^2,size(x,1),1);



end

