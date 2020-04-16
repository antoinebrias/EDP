function [mdlstruct] = init_mdl_struct(n_dim,n_lags)
% INIT_MDL_STRUCT Model structure default  initialization.
%   [MDLSTRUCT] = INIT_MDL_STRUCT(N_DIM) sets a model structure with N_DIM species 
%
%
%   [MDLSTRUCT] = INIT_MDL_STRUCT(N_DIM,N_LAGS) sets a model structure with
%   N_DIM species and N_LAGS lags for each species.
%
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


