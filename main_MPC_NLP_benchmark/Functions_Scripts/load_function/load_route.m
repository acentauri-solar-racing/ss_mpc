%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 01.10.2023
%% Description 
% This function initializes all the variables needed for the mpc simulation
% loop, run the actual simulation and stores the results 

% INPUT: 
% "s_step": discretization step of mpc

% OUTPUT : 
% "route": struct with inclination (alpha), max speed (max_v) and
% cumulative distance (s_x)
%%
function route = load_route(s_step)
%addpath('..\..\ss_offline_data\route\BWSC');
    s_x = 0:s_step:3000000;
    alpha = readtable('route_preprocessed.csv').inclinationSmooth;
    maxSpeed = readtable('route_preprocessed.csv').maxSpeed;
    cumDistance = readtable('route_preprocessed.csv').cumDistance;
    route.alpha = interp1(cumDistance, alpha, s_x);
    route.max_v = interp1(cumDistance, maxSpeed, s_x);
    route.cumDistance = s_x;
end