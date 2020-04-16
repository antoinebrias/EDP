%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Current  value function
% (we look for the control optimizing this function)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ res, resUnweighted ] = fun_optim_td(support_states,control,value_function,model,optstruct,mdlstruct)

discount_factor = optstruct.discount_factor;
reward = optstruct.reward;
weights = optstruct.weights;

n_lags = mdlstruct.n_lags;
ind_available_var = mdlstruct.ind_available_var;
ind_current_var = mdlstruct.ind_current_var;
n_dim=mdlstruct.n_dim;

SS= repelem(support_states,size(control,1),1);
currentU=repmat(control,size(support_states,1),1);
currentX = SS;

[mu,~]=model(currentX,currentU,1); % posterior mean
nextX = mu;%.*currentU;%.*currentU./(currentU+eps);

nextX(nextX<0)=0;

if isempty(value_function)
    [res,resUnweighted]=reward(currentX(:,ind_available_var)*0,currentU*0,weights);

else
  
    [resTotTmp,resUnweighted]=reward(currentX(:,ind_available_var),currentU,weights);
  
    [muV,~]=value_function.gp_model(nextX,1);
    
    resV = muV*weights'; resVI = muV;
    
    res = resTotTmp+discount_factor*resV;
    resUnweighted = resUnweighted+discount_factor*resVI;
    
    
end
%Since we solve the min pb we need the opposite values
res = -res;
resUnweighted = -resUnweighted;
end