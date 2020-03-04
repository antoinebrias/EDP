% GP regression vizualisation

function []=gp_show(mdlstruct)
% name = mdlstruct.name; 

x1l = linspace(0,3,20);
x2l = linspace(0,3,20);
[x1,x2] = meshgrid(x1l,x2l);
M = [x1(:) x2(:)];

 [mu_gp,v_gp] = mdlstruct.gp_model(M,M*0,1);
[mu,v] = mdlstruct.model(M,M*0,1);

figure;surf(x1,x2,reshape(mu(:,1),20,20))
hold on;surf(x1,x2,reshape(mu_gp(:,1),20,20))

plot3(mdlstruct.data(1:end-1,1).*(1-mdlstruct.control_data(1:end-1,1)),mdlstruct.data(1:end-1,2).*(1-mdlstruct.control_data(1:end-1,2)),mdlstruct.data(2:end,1),'rx')



end