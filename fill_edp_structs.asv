% Complete the structures by adding default options if needed, and prepare
% additionnal parameters required for the computation

function [mdlstruct,optstruct,dpstruct]=fill_edp_structs(mdlstruct,optstruct,dpstruct)
%% Model Structure
if isfield(dpstruct,'data')
mdlstruct.n_dim = size(mdlstruct.data,2);
else
    if isfield(dpstruct,'name')
        mdlstruct.n_dim = length(mdlstruct.name);
    end
end






%% Optimal control Structure








%% Dynamic programming Structure

if ~isfield(dpstruct,'method')
dpstruct.method = 'td';
end


if ~isfield(dpstruct,'n_support_states')
% number of support states
dpstruct.n_support_states=50;
end

if dpstruct.method == 'td' & ~isfield(dpstruct,'generating_support_states_method')
% How to generate the support states used here
% 'random' -> randomely generated 
% 'grid'   -> create a grid of support states (number of points on each dimension equal to dpstruct.nb_support_states^(1/n) (round to lower))
dpstruct.generating_support_states_method = 'random';
dpstruct = generate_support_states(mdlstruct,dpstruct);

end

if ~isfield(dpstruct,'lambda')
% TD lambda parameter
dpstruct.lambda=0.5;
end

if ~isfield(dpstruct,'n_out_max')
% outer loop number of iterations
dpstruct.n_out_max = 20;
end

if ~isfield(dpstruct,'n_in_max')
% inner loop number of iterations
dpstruct.n_in_max = 5;
end

if ~isfield(dpstruct,'debug')
% debug mode with additionnal display during the computation
dpstruct.debug = 0;
end




end