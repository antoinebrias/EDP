%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data structure default  initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mdlstruct] = init_mdl_struct(n_dim,n_lags)

mdlstruct.n_dim = n_dim;
mdlstruct.tau = 1;
if nargin >1
    mdlstruct.n_lags=n_lags;
else
    mdlstruct.n_lags=ones(1,n_dim);
end


for i=1:n_dim
    mdlstruct.name{i} =['x_' num2str(i)];
end

end


