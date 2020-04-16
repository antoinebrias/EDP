function [x,u] = generate_data(model,xinit,u,T)
% GENERATE_DATA   Data generation following a given model.
%    [X,U] = GENERATE_DATA(MODEL,XINIT,U,T) generates data according the
%    MODEL, starting from XINIT, with the list of control U (!!! need to add more clever rules), until the time
%    horizon T.

x(1,:) = xinit;


for t = 2:T
    x(t,:)=model(x(t-1,:),u(t-1,:),0);
end


end


