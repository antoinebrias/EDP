%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!! DEPRECATED !!!
% This is a template main file for EDP 1.0.
% The different structures have to be initialized (see edp_structs_template.m)
%
% The 3 data structures are:
%  - mdlstruct  (time-series data)
%  - optstruct  (optimal control problem )
%  - dpstruct   (temporal difference learning)
%
% This script fits a GP regression based on mdlstruct.
% Then, the GP regression is used in the DP algorithm (whose parameters are in dpstruct)
% to solve the optimal control (whose parameters are in optstruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  GP REGRESSION  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this step can be skipped if a model already exists, otherwise a GP is
% fitted to the data
% fills the field mdlstruct.model with the fitted gp
% fills the field mdlstruct.model_stats with statistics computed during the
% GP regression

mdlstruct = fit_gp(mdlstruct);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  APPROX DYNAMIC PROGRAMMING USING  %%
%%  TEMPORAL DIFFERENCE LEARNING  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




