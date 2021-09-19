% Run this script to check/set dependencies in the path
% Third party tools:
%   - manopt
%	--- --- ---


% CVX: modelling tool for convex optimization
if ~exist('dependences/manopt')
  warning('Download manopt and add main folder to path')
else
    % You should have cvx folder inside relative_pose_essential
    run dependences/manopt/importmanopt.m
end


% Add common functions
path_common = genpath(pwd);
rmpath(path_common);
addpath(path_common);


disp('SETUP: All dependencies added correctly')
