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
