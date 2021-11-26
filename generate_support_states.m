function [dpstruct] = generate_support_states(mdlstruct,dpstruct)
% GENERATE_SUPPORT_STATES   Data generation following a given model.
%    [DPSTRUCT] = GENERATE_SUPPORT_STATES(MDLSTRUCT,DPSTRUCT) generates
%    list of support states according the data in MDLSTRUCT and DPSTRUCT
%
%    DPSTRUCT.N_SUPPORT_STATES is the number of support states
%
%    DPSTRUCT.GENERATING_SUPPORT_STATES_METHOD is the method used to
%    generate the list (!!! for now only 'random' is working)



n_lags = mdlstruct.n_lags;
ind_available_var = mdlstruct.ind_available_var;


if ~isfield(dpstruct,'n_support_states')
    % number of support states
    dpstruct.n_support_states=50;
end


if isfield(dpstruct,'n_support_states')
    switch dpstruct.generating_support_states_method
        case 'grid'
            %% to do
            if mdlstruct.n_lags==1
                dpstruct.support_states = linspace( min(mdlstruct.data(:,ind_available_var(1))), max(mdlstruct.data(:,ind_available_var(1))),dpstruct.n_support_states)';
            else
                
                support_states_tmp = [];
                k=1;
                for i=1:mdlstruct.n_dim
                    maxdata = repelem(max(mdlstruct.data(:,ind_available_var(i))),1,n_lags(ind_available_var(i)));
                    mindata = repelem(min(mdlstruct.data(:,ind_available_var(i))),1,n_lags(ind_available_var(i)));
                    support_states_tmp = [support_states_tmp rand(dpstruct.n_support_states,mdlstruct.n_lags(ind_available_var(i)) ).*(maxdata-mindata)+mindata];
                    
                    for j=1:size(mindata,2)
                        current_dim_grid = linspace(mindata(j),maxdata(j),floor(sqrt(dpstruct.n_support_states)));
                        grids{k} = current_dim_grid;
                        x_grids{k}={};
                        k=k+1;
                    end
                end
                
                [x_grids{:}] = meshgrid(grids{:});
                M = [x_grids{1}(:) x_grids{2}(:)];
                
                
                dpstruct.support_states =M;
                
                dpstruct.n_support_states=floor(sqrt(dpstruct.n_support_states))^2;
                
                
            end
            
        case 'random'
            support_states_tmp = [];
            for i=1:mdlstruct.n_dim
                maxdata = repelem(max(mdlstruct.data(:,ind_available_var(i))),1,n_lags(ind_available_var(i)));
                mindata = repelem(min(mdlstruct.data(:,ind_available_var(i))),1,n_lags(ind_available_var(i)));
                support_states_tmp = [support_states_tmp rand(dpstruct.n_support_states,mdlstruct.n_lags(ind_available_var(i)) ).*(maxdata-mindata)+mindata];
                
            end
            dpstruct.support_states =support_states_tmp;
    end
end