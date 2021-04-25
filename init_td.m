%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Temporal difference learning structure default  initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dpstruct] = init_td(mdlstruct)
% INIT_TD Temporal difference learning structure default  initialization.
%   [DPSTRUCT] = INIT_TD(MDLSTRUCT) sets a default temporal difference
%   learning (approximate dynamic programming) structure. MDLSTRUCT is an
%   input in order to generate the support states list.
%

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
dpstruct.lambda=0.9;

% outer loop number of iterations
dpstruct.n_out_max = 50;

% inner loop number of iterations
dpstruct.n_in_max = 1;

% debug mode with additionnal display during the computation
dpstruct.debug = 0;

%embedded dimension information
dpstruct.ind_available_var = mdlstruct.ind_available_var;


end







