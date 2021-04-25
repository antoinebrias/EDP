% Time-series upload
clear all
close all
rng('default');
rng(1);


% data = load('C:\Users\Renaud\Documents\MATLAB\edp2\ricker1.mat');
 data = load('C:\Users\Renaud\Documents\MATLAB\edp2\sailalorda1.mat');




x = data.xt(1,:)';
u = data.u';


% Time-series length
T=size(x,1);


% default EDM structure
[mdlstruct] = init_mdl_struct(1); % the input is the number of available species data

% data generation
mdlstruct.data = x;
mdlstruct.control_data = u ;

%different kinds of control exist, we need to specify which one we use here
mdlstruct.control_type =  'rate';
mdlstruct.control_param =  [1 ]; % both  dimensions are fully controlled

% data display
data_show(mdlstruct)


% we specify the number of lags used here
% here, we specify the number of time lags used for each variables
% 0 means that the variable is not used in the prediction
mdlstruct.n_lags = [2];

% GP regression
disp('*****************************************************')
disp('GP fitting')
disp('*****************************************************')
% this step can be skipped if a model already exists, otherwise a GP is
% fitted to the data
% fills the field mdlstruct.model with the fitted gp
% fills the field mdlstruct.model_stats with statistics computed during the
% GP regression

% default gp structure
mdlstruct.gp = init_gp();

mdlstruct.gp.is_log = 0;
mdlstruct.gp.cond_0=1;

% we now fit the gp with the data in mdlstruct
mdlstruct = fit_gp(mdlstruct);

% gp display
gp_show(mdlstruct)


% Optimal control problem
disp('*****************************************************')
disp('  Optimal control problem structure creation')
disp('*****************************************************')

% Discount factor
optstruct.discount_factor=0.95;

% Objective weights
optstruct.weights = [1];


% Reward function
% The instantaneous reward function is a handle function.
% The input are the current variables, the controls and the weight for each
% dimension
% The output are the total weighted reward, the unweighted reward on each
% dimension, and the weighted reward on each dimension

%>>> Example 1 <<<<
% use of an existing function
optstruct.reward = @(x,u,w) reward_harvest_all(x,u,w);

% Temporal difference learning
disp('*****************************************************')
disp('EDP using Temporal difference learning')
disp('*****************************************************')
% default  td structure
dpstruct=init_td(mdlstruct);

%  run the td algorithm, using the gp regression
 [dpstruct] = td_learning(mdlstruct,optstruct,dpstruct,'gp');

 % control maps
 disp('Control maps given by the EDP algorithm')
td_show(mdlstruct,optstruct,dpstruct)

% controlled trajectories
t_max = 50; % time length
n_traj = 10; % number of trajectories
start_x = 'last'; % starting point
is_disp = 1; % display the trajectories result
 disp('Controlled trajectories following the EDP policy')
traj = sim_trajectories(mdlstruct,optstruct,dpstruct,t_max,n_traj,start_x,is_disp);

% Pareto front
disp('*****************************************************')
disp('Pareto front')
disp('*****************************************************')
n_pareto = 10;
pareto = pareto_front(mdlstruct,optstruct,dpstruct,n_pareto,'gp');

% % control map for different Pareto weights
% td_show(mdlstruct,pareto.optstruct{1},pareto.dpstruct{1})
% td_show(mdlstruct,pareto.optstruct{end},pareto.dpstruct{end})






function y = competition(x,u,r)
    r2 = 2/3*r; m1 = 2; m2 = 1; b = 4;
    s = x.*(1-u);
    y = [r r2].*exp(normrnd(-0.01/2,0.1)).*s.*[1+(s(:,1)+m1*s(:,2)).^b  1+(s(:,2)+m2*s(:,1)).^b].^-1;
end

function y = migration(x,u,r)
    r2 = 2/3*r; m1 = 2; m2 = 1; b = 4;
    s = x.*(1-u);
    y = [r r2].*s.*[1+(s(:,1)+m1*s(:,2)).^b  1+(s(:,2)+m2*s(:,1)).^b].^-1;
end

function y = maternal(x,u,r)
    g = 10; c=5; 
    s = x.*(1-u);
    y = [r r2].*s.*[1+(s(:,1)+m1*s(:,2)).^b  1+(s(:,2)+m2*s(:,1)).^b].^-1;
end