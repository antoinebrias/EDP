%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Data generation following a given model
%
% Input are: the model used (model(x,u,is_det)), the initial starting point xinit, the
% control u (if available), the length of the time series, the available
% variables
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,u] = generate_data(model,xinit,u,T)
x(1,:) = xinit;


for t = 2:T
    x(t,:)=model(x(t-1,:),u(t-1,:),0);
end


end


