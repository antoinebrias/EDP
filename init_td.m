%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Temporal difference learning structure default  initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dpstruct] = init_td(mdlstruct)
dpstruct.method = 'td';

%% Backward dynamic programming parameters


%% temporal difference learning parameters

% number of support states
dpstruct.n_support_states=20;

% How to generate the support states used here
% 'random' -> randomely generated 
% 'grid'   -> create a grid of support states (number of points on each dimension equal to dpstruct.nb_support_states^(1/n) (round to lower))
dpstruct.generating_support_states_method = 'random';
if nargin==1
dpstruct = generate_support_states(mdlstruct,dpstruct);
end
% TD lambda parameter
dpstruct.lambda=0.5;

% outer loop number of iterations
dpstruct.n_out_max = 10;

% inner loop number of iterations
dpstruct.n_in_max = 2;

% debug mode with additionnal display during the computation
dpstruct.debug = 0;


end







