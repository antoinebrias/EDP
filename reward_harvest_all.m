%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reward function: harvest all
% Reward on every species catch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [total_reward,unweighted_reward,weighted_reward] = reward_harvest_all(x,u,weights)
unweighted_reward=[x.*u ]; % reward on the harvest of x_1 and escapement of x_2
weighted_reward = unweighted_reward.*weights;   % the objectives are weighted
total_reward = sum(weighted_reward,2);  % total reward

end


