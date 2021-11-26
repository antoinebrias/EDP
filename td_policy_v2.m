function [opt_control,unweighted_opt_value,weighted_opt_value]=td_policy_v2(x,optstruct,dpstruct,mdlstruct)
% TD_POLICY  finds the optimal control given the TD value function and the
% state of the system
%    [OPT_CONTROL,UNWEIGHTED_OPT_VALUE,WEIGHTED_OPT_VALUE]=TD_POLICY(X,OPTSTRUCT,DPSTRUCT,MDLSTRUCT)
%    returns the optimal control according the EDP policy for the system in
%    state X.


%% One-step ahead DP problem
% solve the one-step ahead DP equation by looking for the control
% maximizing fun_optim_td
fOptTD = @(control) fun_optim_td_v2 (x,control,dpstruct.value_function,dpstruct.model,optstruct,mdlstruct);

ftd = [];

[ftd,~]=fOptTD(dpstruct.control_grid);


[~,indCurrentU]=min(reshape(ftd,size(dpstruct.control_grid,1),size(x,1)));

% optimal control
opt_control=dpstruct.control_grid(indCurrentU,:);

[muV,~]=eval_value_function(dpstruct.value_function,x);
% [muV,~]=dpstruct.value_function.gp_model(x,1);
weighted_opt_value = muV*optstruct.weights';
unweighted_opt_value = muV;


end