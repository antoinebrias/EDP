%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Support states list generation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dpstruct] = generate_support_states(mdlstruct,dpstruct)
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
        case 'random'
            support_states_tmp = [];
            for i=1:mdlstruct.n_dim
            maxdata = repelem(max(mdlstruct.data(:,ind_available_var(i))),1,n_lags(ind_available_var(i)));
            mindata = repelem(min(mdlstruct.data(:,ind_available_var(i))),1,n_lags(ind_available_var(i)));
            support_states_tmp = [support_states_tmp rand(dpstruct.n_support_states,mdlstruct.n_lags(ind_available_var(i)) ).*(maxdata-mindata)-mindata];
            
            end
             dpstruct.support_states =support_states_tmp;
end
end
