function [total_reward,unweighted_reward,weighted_reward] = reward_harvest_all(x,u,weights)
% REWARD_HARVEST_ALL  Reward function: reward on every species catched.
%   [TOTAL_REWARD,UNWEIGHTED_REWARD,WEIGHTED_REWARD] = REWARD_HARVEST_ALL(X,U,WEIGHTS)
%   Return the reward, applying the control U on the current state X, with
%   a weight WEIGHTS on each variable.
%
%   TOTAL_REWARD is the sum of every species rewards. UNWEIGHTED_REWARD is
%   a vector containing each individual reward.  WEIGHTED_REWARD is
%   a vector containing each weighted individual reward.
%

unweighted_reward=[x.*u ]; % reward on the harvest of x
weighted_reward = unweighted_reward.*weights;   % the objectives are weighted
total_reward = sum(weighted_reward,2);  % total reward

end


