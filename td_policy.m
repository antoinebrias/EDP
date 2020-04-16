%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% find the optimal control given the TD value function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [opt_control,unweighted_opt_value,weighted_opt_value]=td_policy(x,optstruct,dpstruct,mdlstruct)

%% One-step ahead DP problem
% solve the one-step ahead DP equation by looking for the control
% maximizing fun_optim_td
fOptTD = @(control) fun_optim_td (x,control,dpstruct.value_function,dpstruct.model,optstruct,mdlstruct);

ftd = [];

[ftd,~]=fOptTD(dpstruct.control_grid);


[~,indCurrentU]=min(reshape(ftd,size(dpstruct.control_grid,1),size(x,1)));

% optimal control
opt_control=dpstruct.control_grid(indCurrentU,:);

[muV,~]=dpstruct.value_function.gp_model(x,1);
weighted_opt_value = muV*optstruct.weights';
unweighted_opt_value = muV;


end