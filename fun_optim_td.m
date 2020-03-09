%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Current  value function
% (we look for the control optimizing this function)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ res, resUnweighted ] = fun_optim_td(support_states,control,value_function,model,optstruct,dpstruct)

discount_factor = optstruct.discount_factor;
reward = optstruct.reward;
weights = optstruct.weights;



SS= repelem(support_states,size(control,1),1);
currentU=repmat(control,size(support_states,1),1);
currentX = SS;

[mu,~]=model(currentX,currentU,1); % posterior mean
nextX = mu;%.*currentU;%.*currentU./(currentU+eps);

nextX(nextX<0)=0;

if isempty(value_function)
    [res,resUnweighted]=reward(currentX*0,currentU*0,weights);

else
  
    [resTotTmp,resUnweighted]=reward(currentX,currentU,weights);
  
    [muV,~]=value_function.gp_model(nextX,1);
    
    resV = muV*weights'; resVI = muV;
    
    res = resTotTmp+discount_factor*resV;
    resUnweighted = resUnweighted+discount_factor*resVI;
    
    
end
%Since we solve the min pb we need the opposite values
res = -res;
resUnweighted = -resUnweighted;
end