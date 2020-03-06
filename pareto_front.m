%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pareto front c omputation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pareto = pareto_front(mdlstruct,optstruct,dpstruct,n_pareto)


if mdlstruct.n_dim == 2
    rho  = linspace(0,1,n_pareto);
    Rho=[rho' 1-rho'];
    
    
        optstruct.weights =[1 0];
    [dpstruct10] = td_learning(mdlstruct,optstruct,dpstruct,'gp');
    
    traj10 = sim_trajectories(mdlstruct,optstruct,dpstruct,50,10,'last');
    
    
        optstruct.weights = [0 1];
    [dpstruct01] = td_learning(mdlstruct,optstruct,dpstruct,'gp');
       traj01 = sim_trajectories(mdlstruct,optstruct,dpstruct,50,10,'last');
    
    
    for i=1:n_pareto
    optstruct.weights = Rho(i,:)./[unique(traj10.th_value) unique(traj01.th_value)];
    [dpstruct] = td_learning(mdlstruct,optstruct,dpstruct,'gp')
     traj{i} = sim_trajectories(mdlstruct,optstruct,dpstruct,50,10,'last');
    end
    
    pareto.traj = traj;
    
end
% 
% % pareto weights
% rho1=linspace(0,1,n_pareto)';rho2=linspace(0,1,n_pareto)';
% [Rho1,Rho2]=meshgrid(rho1,rho2);
% Mrho = [Rho1(:) Rho2(:)];
% Mrho(sum(Mrho,2)>1,:)=[];
% Mrho(:,3) = 1-sum(Mrho,2);
% nrho = size(Mrho,1);
% 
% 
% 
% 
% % temporal difference learning
% 
% % default  td structure
% dpstruct=init_td(mdlstruct);
% 
% %  run the td algorithm, using the gp regression
%  [dpstruct] = td_learning(mdlstruct,optstruct,dpstruct,'mdl');
% 
%  % control maps
% td_show(mdlstruct,optstruct,dpstruct)
% 
% % controlled trajectories
% t_max = 50;
% n_traj = 10;
% start_x = 'last';
% traj = sim_trajectories(mdlstruct,optstruct,dpstruct,t_max,n_traj,start_x);
% 
% % Pareto front
% n_pareto = 10;
% pareto = pareto_front(mdlstruct,optstruct,dpstruct,n_pareto);