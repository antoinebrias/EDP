%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate controlled trajectories
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function traj = sim_trajectories(mdlstruct,optstruct,dpstruct,t_max,n_traj,start_x,is_display)

if nargin <7
    is_display = 0;
end

n_dim = mdlstruct.n_dim;

if ischar(start_x)
    switch start_x
        case 'last'
            x_init = mdlstruct.data(end,:);
        case 'mean'
            x_init = mean(mdlstruct.data);
    end
else
    x_init = start_x;
end

X{1} = repmat(x_init, n_traj,1);
tot_value = 0;
tot_unweighted_value=zeros(1,n_dim);
for t=2:t_max
    
    [opt_control_tmp,unweighted_opt_value,weighted_opt_value]=td_policy(X{t-1},optstruct,dpstruct);
    
    if  t==2
        traj.th_value = weighted_opt_value;
        traj.th_unweighted_value = unweighted_opt_value;
    end
    
    
    opt_control{t-1} = opt_control_tmp;
    if isfield(mdlstruct,'model')
        X{t} = mdlstruct.model(X{t-1},opt_control{t-1},0);
    else
        if isfield(mdlstruct,'gp_model')
            X{t} =mdlstruct.gp_model(X{t-1},opt_control{t-1},0);
        end
    end
    
    
    [instant_value_tmp{t-1} unweighted_instant_value_tmp{t-1}] = optstruct.reward(X{t-1},opt_control{t-1},optstruct.weights);
    
    tot_unweighted_value =  tot_unweighted_value+unweighted_instant_value_tmp{t-1}*optstruct.discount_factor^(t-2);
    tot_value = tot_value+instant_value_tmp{t-1}*optstruct.discount_factor^(t-2);
    instant_cum_value(t,1:n_traj) = tot_value;
end


for t=1:t_max-1
    data(t,1:n_traj,1:n_dim) = X{t};
    control(t,1:n_traj,1:n_dim) =opt_control{t};
end



traj.data = data;
traj.control = control;
traj.value = tot_value;
traj.instant_cum_value = instant_cum_value;
traj.unweighted_value =  tot_unweighted_value;
traj.weights =optstruct.weights;


name = mdlstruct.name;



%% display
if is_display
    
    mean_x = squeeze(mean(traj.data));
    mean_catch = squeeze(mean(traj.data.*traj.control));
    mean_v = squeeze(mean(traj.instant_cum_value'));
    figure
    title('Mean trajectory')
    hold on
    [pplot,nplot]=numSubplots(n_dim+1);
    for i=1:n_dim
        subplot(pplot(1),pplot(2),i);
        hold on
        plot(mean_x(:,i),'b','Linewidth',2);
        
        switch mdlstruct.control_type
            case 'rate'
                plot(mean_catch(:,i),'c','Linewidth',2);
                leg_control = [name{i} ' catch'];
            case 'single'
                % to do
            case 'global'
                % to do
                
        end
        
        legend(name{i},leg_control);
        xlabel('Time')
        ylabel(name{i})
        
        grid on;  box on
    end
    
    subplot(pplot(1),pplot(2),n_dim+1);
    hold on
    plot(mean_v,'g','Linewidth',2);
    legend('cumulative reward');
    xlabel('Time')
    ylabel('Reward')
    grid on;  box on
    
end



end