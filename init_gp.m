function [gp] = init_gp()
% INIT_GP  Gaussian Process structure default  initialization.
%   [GP] = INIT_GP() sets a GP regression up with default options 


%% GP options
% by default the regression is Y=GP(X) but if we want Y=X*exp(GP(X)), this
% parameter should be set to 1
gp.is_log=0;

%  set to 1 to put a condition on the origin
gp.cond_0=0;

% set to [j 0] to put a condition when the variable j is 0
gp.cond_j_0=[0 0];

% set to [j 1] to put a condition when the variable j is 1
gp.cond_j_1=[0 0];



end







