%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Support states list generation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dpstruct] = generate_support_states(mdlstruct,dpstruct)
n_lags_max = mdlstruct.n_lags_max;

if ~isfield(dpstruct,'n_support_states')
    % number of support states
    dpstruct.n_support_states=50;
end


if isfield(dpstruct,'n_support_states')
    switch dpstruct.generating_support_states_method
        case 'grid'
            %% to do
        case 'random'
            maxdata = repelem(max(mdlstruct.data),1,n_lags_max);
            mindata = repelem(min(mdlstruct.data),1,n_lags_max);
            dpstruct.support_states = rand(dpstruct.n_support_states,mdlstruct.n_lags_max.*mdlstruct.n_dim ).*(maxdata-mindata)-mindata;
            
            
    end
else
    % default case
    maxdata = repelem(max(mdlstruct.data),1,n_lags_max);
    mindata = repelem(min(mdlstruct.data),1,n_lags_max);
    dpstruct.support_states = rand(dpstruct.n_support_states,mdlstruct.n_lags_max.*mdlstruct.n_dim ).*(maxdata-mindata)-mindata;
    dpstruct.generating_support_states_method ='random';
end
end
