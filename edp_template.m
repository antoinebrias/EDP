%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a template inputfile for EDP 1.0.
% The data are given in 3 mandatory data structures, with various data fields.
%
% The 3 data structures are:
%  - mdldata  (time-series data)
%  - optdata  (optimal control problem )
%  - dpdata   (temporal difference learning)
%
% The number of data fields in these structures is variable and some of them may not have to be specified.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'MDLSTRUCT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.
% Use no spaces in names.
% probdata.name =  { 'name1' 'name2' ... } or { 'name1' 'name2' ... }'

%>>> Example 1 <<<<
mdlstruct.name =  { 'X1'
                   'X2'
                   'X3'
                   'X4'};

%>>> Example 2 <<<<
mdlstruct.name =  { 'catch'
                   'cpue'};

% Type of model used: 
% 'gp'  -> GP regression
% 'fct' -> handle function
mdlstruct.mdl_type='gp';
mdlstruct.mdl_type='fct';

% example of function
mdlstruct.fct = @(x,u,is_det) ricker2d(x,u,is_det,param) ;



% Control variables
% Specify which variable can be controlled

%>>> Example 1 <<<<
mdlstruct.control =  { 'X1'
                   'X4'};

% %>>> Example 2 <<<<
% mdlstruct.control =  [1 0 0 1]; % the first and fourth variables are controlled



% Control type
% Specify how the variables are controlled
% control_param, is an additionnal parameter for the behaviour of the
% control
% 'rate"  ->  independent fraction of the controlled variables is removed (between 0 and 1)
% 'single" -> one control for several variables (bycatch for example)
%             control_param contains the bycatch rate for each controlled
%             variable
% 'global" ->   TO DO
%             control_param contains the min and max possible control on
%             each variable


%>>> Example 1 <<<<
mdlstruct.control_type =  'rate';

%>>> Example 2 <<<<
mdlstruct.control_type =  'single';
mdlstruct.control_param =  [1 0.2 ]; % here, the escapement is x_1.*(1-u) and  x_2.*(1-0.2*u) 

%>>> Example 3 <<<<
mdlstruct.control_type =  'global';
mdlstruct.control_param =  [0 10; 0 5]'; % optionnal, x_1 can be controlled between 0 and 10, x_2 between 0 and 5





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'OPTSTRUCT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discount factor
optstruct.discount_factor=0.95;

% Objective weights
optstruct.weights = [0.5 0.5];


% Reward function
% The instantaneous reward function is a handle function.
% The input are the current variables, the controls and the weight for each
% dimension
% The output are the total weighted reward, the unweighted reward on each
% dimension, and the weighted reward on each dimension

%>>> Example 1 <<<<
% use of an existing function
optstruct.reward = @(x,u) harvest_all(x,u,optstruct.weights)


%>>> Example 2 <<<<
% local function at the end of this script
optstruct.reward = @(x,u) local_reward(x,u,optstruct.weights);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'DPSTRUCT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Two methods can be used, the backward dynamic programming (for small
% problem) and the temporal difference learning
% 'bdp' ->  backward dynamic programming method
% 'td'  ->  temporal difference learning

dpstruct.method = 'td';

%% Backward dynamic programming parameters


%% temporal difference learning parameters

% number of support states
dpstruct.nb_support_states=50;

% How to generate the support states used here
% 'random' -> randomely generated 
% 'grid'   -> create a grid of support states (number of points on each dimension equal to dpstruct.nb_support_states^(1/n) (round to lower))
dpstruct.generating_support_states_method = 'random';


% TD lambda parameter
dpstruct.lambda=0.5;

% outer loop number of iterations
dpstruct.nOutMax = 20;

% inner loop number of iterations
dpstruct.nInMax = 5;

% debug mode with additionnal display during the computation
dpstruct.debug = 0;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS COMPLETION  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function fills the structures to be ready to use in the edp
% algorithm. It contains all the default options.
[mdlstruct,optstruct,dpstruct]=fill_edp_structs(mdlstruct,optstruct,dpstruct);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LOCAL FUNCTIONS  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example of reward function
function [total_reward,unweighted_reward,weighted_reward]=local_reward(x,u,weights)
unweighted_reward=[x(:,1).*u(:,1) x(:,2).*(1-u(:,2)) ]; % reward on the harvest of x_1 and escapement of x_2
weighted_reward = unweighted_reward.*weights;   % the objectives are weighted
total_reward = sum(resWeighted,2);  % total reward
end


