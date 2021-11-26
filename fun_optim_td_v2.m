function [res, resUnweighted] = fun_optim_td_v2(support_states,control,value_function,model,optstruct,mdlstruct)
% FUN_OPTIM_TD   Value function of the one-step ahead Dynamic Programming
% equation, which has to be optimized according the control.
%    [RES, RESUNWEIGHTED] = FUN_OPTIM_TD(SUPPORT_STATES,CONTROL,VALUE_FUNCTION,MODEL,OPTSTRUCT,MDLSTRUCT) fits
%    We are looking for the control optimizing this function. 
%    For each state in the SUPPORT_STATES list, we apply every controls and
%    return the result of the one-step ahead Dynamic Programming
%    equation. RES contains the weighted sum  over all objectives. RESUNWEIGHTED
%    contains the result for each objective.

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



update_u=[];
for i=1:n_dim
    utmp = [currentU(:,i) zeros(size(currentU,1),n_lags(ind_available_var(i))-1)];
    update_u=[update_u utmp];
end
currentU = update_u;



[mu,~]=model(currentX,currentX.*0,1); % posterior mean
% nextX = mu;%.*currentU;%.*currentU./(currentU+eps);

% nextX(nextX<0)=0;

 nextX=[];
        for i=1:n_dim
            if n_lags(i)==1
                % one lag: nextX = mu
                 xtmp = [mu(:,i).*(1-currentU(:,ind_available_var(i)))];
%                  xtmp = [mu(:,i) currentX(:,ind_available_var(i)).*(1-currentU(:,ind_available_var(i)))];
            else
                 if n_lags(i)==2
                     % 2 lags: nextX =[mu X(1-U)]
%                  xtmp = [mu(:,i) currentX(:,ind_available_var(i)).*(1-currentU(:,ind_available_var(i)))];
                 
                 xtmp = [mu(:,i).*(1-currentU(:,ind_available_var(i))) currentX(:,ind_available_var(i))];
                 
                 
%                 xtmp = [mu(:,i)  currentX(:,ind_available_var(i)).*(1-currentU(:,ind_available_var(i))) currentX(:,ind_available_var(i)+1:ind_available_var(i)+n_lags(i)-1)];
                 else
                     % number of lags over 2: nextX =[mu X(1-U) S ....]
%                      xtmp = [mu(:,i)  currentX(:,ind_available_var(i)).*(1-currentU(:,ind_available_var(i))) currentX(:,ind_available_var(i)+1:ind_available_var(i)+n_lags(i)-1)];
                    
                      xtmp = [mu(:,i).*(1-currentU(:,ind_available_var(i)))  currentX(:,ind_available_var(i)) currentX(:,ind_available_var(i)+1:ind_available_var(i)+n_lags(i)-1)];
                      
                     
                     
                     
                 end
            end
            nextX=[nextX xtmp];
        end
        
%         nextX(nextX<0)=0;

if isempty(value_function)
    [res,resUnweighted]=reward(nextX(:,ind_available_var)*0,currentU*0,weights);

else
  
    [resTotTmp,resUnweighted]=reward(mu(:,ind_available_var),currentU(:,ind_available_var),weights);
  
%     [muV,~]=value_function.gp_model(nextX,1);
    [muV,~]=eval_value_function(value_function,nextX);
    
    resV = muV*weights'; resVI = muV;
    
    res = resTotTmp+discount_factor*resV;
    resUnweighted = resUnweighted+discount_factor*resVI;
    
    
end
%Since we solve the min pb we need the opposite values
res = -res;
resUnweighted = -resUnweighted;
end