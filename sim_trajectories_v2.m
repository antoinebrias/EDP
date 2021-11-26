function traj = sim_trajectories_v2(mdlstruct,optstruct,dpstruct,t_max,n_traj,start_x,is_display)
%  V2 CHANGES !!!!!!
% SIM_TRAJECTORIES  V2  Simulate controled trajectories according EDP results
%   SIM_TRAJECTORIES(MDLSTRUCT,OPTSTRUCT,DPSTRUCT,T_MAX,N_TRAJ,START_X,IS_DISPLAY)
%   simulates N_TRAJ starting from START_X, until time horizon T_MAX. At
%   each step, the optimal control is computed thanks to the EDP policy.
%   START_X may be a string: 'mean' for the mean point in the dataset,
%   'last' for the last point in the dataset, or a vector.
%
%   IS_DISPLAY indicates if a figure containing the mean trajectory is produced.

n_lags = mdlstruct.n_lags;
ind_available_var = mdlstruct.ind_available_var;
ind_current_var = mdlstruct.ind_current_var;
n_dim=mdlstruct.n_dim;

if nargin <7
    is_display = 0;
end

n_dim = mdlstruct.n_dim;
n_real_dim=size(mdlstruct.data,2);

x_init = [];
if ischar(start_x)
    switch start_x
        case 'last'
            if max(n_lags)==1
                x_init = mdlstruct.data(end,:);
            else
                for i= 1:size(mdlstruct.data,2)
                    x_init = [x_init mdlstruct.data(end:-1:end-n_lags(i)+1,:)'];
                end
            end
        case 'mean'
            if max(n_lags)==1
                x_init = mean(mdlstruct.data);
            else
                for i= 1:size(mdlstruct.data,2)
                    x_init = [x_init mean(mdlstruct.data).*ones(1,n_lags(i))];
                end
            end
    end
else
    x_init = start_x;
end

X{1} = repmat(x_init, n_traj,1);
tot_value = 0;
tot_unweighted_value=zeros(1,n_dim);
for t=2:t_max
    t
    [opt_control_tmp,unweighted_opt_value,weighted_opt_value]=td_policy_v2(X{t-1},optstruct,dpstruct,mdlstruct);
    
    if  t==2
        traj.th_value = weighted_opt_value;
        traj.th_unweighted_value = unweighted_opt_value;
    end
    
    
    opt_control{t-1} = opt_control_tmp;
    if isfield(mdlstruct,'model')
        % in this case the control must be set up to 0 for each variables
        % unused in the td learning
        
        tmp_control = zeros(n_traj,n_real_dim);
        tmp_control(:,ind_available_var) =    opt_control{t-1};
        
        %         X{t} = mdlstruct.model(X{t-1},tmp_control,0);
        
        tmp_X = mdlstruct.model(X{t-1},tmp_control,0);
        next_X = [];
        for i= 1:size(mdlstruct.data,2)
            %                  if n_lags(i)>1
            %                      next_X = [next_X tmp_X(:,i) X{t-1}(:,ind_current_var(i):ind_current_var(i)+n_lags(i)-2)];
            %                  else
            %                      next_X = [next_X tmp_X(:,i)];
            %                  end
            
            if n_lags(i)==1
                next_X = [next_X mu(:,i) tmp_X(:,i) X{t-1}(:,ind_available_var(i)).*(1-tmp_control(:,ind_available_var(i)))];
            else
                if n_lags(i)>1
                    next_X = [next_X mu(:,i) tmp_X(:,i)  X{t-1}(:,ind_available_var(i)).*(1-tmp_control(:,ind_available_var(i))) X{t-1}(:,ind_available_var(i)+1:ind_available_var(i)+n_lags(i)-2)];
                end
            end
            
            
            
        end
        X{t} = next_X;
        
        
    else
        if isfield(mdlstruct,'gp_model')
            %             X{t} =mdlstruct.gp_model(X{t-1},opt_control{t-1},0);
            tmp_X =mdlstruct.gp_model(X{t-1},opt_control{t-1},0);
            tmp_X(tmp_X<0) = 0;
            next_X = [];
            tmp_control = opt_control{t-1};
            for i= 1:size(mdlstruct.data,2)
                %                  if n_lags(i)>1
                %                      next_X = [next_X tmp_X(:,i) X{t-1}(:,ind_current_var(i):ind_current_var(i)+n_lags(i)-2)];
                %
                %                  else
                %                      next_X = [next_X tmp_X(:,i)];
                %                  end
                
                if n_lags(i)==1
                    next_X = [next_X  tmp_X(:,i)];
%                      next_X = [next_X  tmp_X(:,i) X{t-1}(:,ind_available_var(i)).*(1-tmp_control(:,ind_available_var(i)))];
                else
                    if n_lags(i)== 2
                         next_X = [next_X  tmp_X(:,i) X{t-1}(:,ind_available_var(i)).*(1-tmp_control(:,ind_available_var(i)))];
                    else
                        % nlag>2
%                            next_X = [next_X  tmp_X(:,i) X{t-1}(:,ind_available_var(i)).*(1-tmp_control(:,ind_available_var(i)))];
                        next_X = [next_X tmp_X(:,i)  X{t-1}(:,ind_available_var(i)).*(1-tmp_control(:,ind_available_var(i))) X{t-1}(:,ind_available_var(i)+1:ind_available_var(i)+n_lags(i)-2)];
                    end
                end
                
                
                
            end
            X{t} = next_X;
            
        end
    end
    
    
    [instant_value_tmp{t-1} unweighted_instant_value_tmp{t-1}] = optstruct.reward(X{t-1}(:,ind_available_var),opt_control{t-1},optstruct.weights);
    
    tot_unweighted_value =  tot_unweighted_value+unweighted_instant_value_tmp{t-1}*optstruct.discount_factor^(t-2);
    tot_value = tot_value+instant_value_tmp{t-1}*optstruct.discount_factor^(t-2);
    instant_cum_value(t,1:n_traj) = tot_value;
end


for t=1:t_max-1
    data(t,1:n_traj,1:n_real_dim) = X{t}(:,ind_current_var);
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
    
    mean_x = squeeze(mean(traj.data,2));
    mean_catch = squeeze(mean(traj.data(:,:,ind_available_var).*traj.control,2));
    if n_traj>1
    mean_v = squeeze(mean(traj.instant_cum_value'));
    else
            mean_v = squeeze(traj.instant_cum_value');
    end
    figure
    title('Mean trajectory')
    hold on
    [pplot,nplot]=numSubplots(n_real_dim+1);
    for i=1:n_real_dim
        subplot(pplot(1),pplot(2),i);
        hold on
        plot(mean_x(:,i),'b','Linewidth',2);
        
        switch mdlstruct.control_type
            case 'rate'
                if n_lags(i) ~= 0
                    plot(mean_catch(:,i),'c','Linewidth',2);
                    leg_control = [name{i} ' catch'];
                end
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
    
    subplot(pplot(1),pplot(2),n_real_dim+1);
    hold on
    plot(mean_v,'g','Linewidth',2);
    %     legend('cumulative reward');
    xlabel('Time')
    ylabel('Reward')
    grid on;  box on
    
end


disp('**********************************************************')
disp(['Theorical value: ' num2str(mean(traj.th_value))])
disp(' ')
disp(['Simulated mean value: ' num2str(mean(traj.value))])
disp(' ')
disp('**********************************************************')

end