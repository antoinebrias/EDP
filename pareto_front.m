function pareto = pareto_front(mdlstruct,optstruct,dpstruct,n_pareto,name_model)
% PARETO_FRONT  Pareto Front computation.
%   PARETO = PARETO_FRONT(MDLSTRUCT,OPTSTRUCT,DPSTRUCT,N_PARETO,NAME_MODEL)
%   Computes the structure PARETO containing the result of EDP for a range
%   of value weightings.
%
% !!! handles 2d case only



if mdlstruct.n_dim == 2
    rho  = linspace(0,1,n_pareto);
    Rho=[rho' 1-rho'];
    
    
    optstruct.weights =[1 0];
    [dpstruct10] = td_learning(mdlstruct,optstruct,dpstruct,name_model);
    
    traj10 = sim_trajectories(mdlstruct,optstruct,dpstruct,50,10,'last');
    
    
    optstruct.weights = [0 1];
    [dpstruct01] = td_learning(mdlstruct,optstruct,dpstruct,name_model);
    traj01 = sim_trajectories(mdlstruct,optstruct,dpstruct,50,10,'last');
    %
    
    for i=1:n_pareto
        optstruct.weights = Rho(i,:)./[unique(traj10.th_value) unique(traj01.th_value)];
        [dpstruct] = td_learning(mdlstruct,optstruct,dpstruct,name_model);
        traj{i} = sim_trajectories(mdlstruct,optstruct,dpstruct,50,10,'last');
        
        
        th_pareto_f(i,:) =  mean(traj{i}.th_unweighted_value);
        pareto_f(i,:) = mean(traj{i}.unweighted_value);
        
        pareto.dpstruct{i} = dpstruct;
        pareto.optstruct{i} = optstruct;
    end
    
    pareto.traj = traj;
    
    
    
    
    figure
    
    
    
    m1=1;%rho'./max(Vrho1);
    m2=1;%(1-rho)'./max(Vrho0);
    
    subplot(1,2,1)
    plot(th_pareto_f(:,1)./m1,th_pareto_f(:,2)./m2,'o','MarkerSize',9,'MarkerFaceColor','b','MarkerEdgeColor','k');hold on;
    
    
    xlabel([mdlstruct.name{1}  ' Long Run Discounted Yield'])
    ylabel({[mdlstruct.name{2} ' Long Run'],'Discounted Yield'})
    title({'Expected Long Term Reward'})
    
    subplot(1,2,2)
    plot(pareto_f(:,1)./m1,pareto_f(:,2)./m2,'o','MarkerSize',9,'MarkerFaceColor','b','MarkerEdgeColor','k');hold on;
    
    
    xlabel([mdlstruct.name{1}  ' Long Run Discounted Yield'])
    ylabel({[mdlstruct.name{2} ' Long Run'],'Discounted Yield'})
    title({'Simulated Long Term Reward'})
    
    
    
    
    pareto.th_pareto_f =  th_pareto_f;
    pareto.pareto_f = pareto_f;
    
end




end
