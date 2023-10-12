%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 01.10.2023
%% Description 
% This function load into the parameters file the route parameters for the
% whole race (from 0 to 3'000 km/3'000'000 m)
% "par.route.incl" = inclination

% INPUT: 
% "par": parameters struct

% OUTPUT : 
% "par": parameters struct with added route parameters
%%
function par = load_route(par)
    s_x = 0:par.s_step:3000000;
    incl = readtable('route_preprocessed.csv').inclinationSmooth;
    maxSpeed = readtable('route_preprocessed.csv').maxSpeed/3.6;
    cumDistance = readtable('route_preprocessed.csv').cumDistance;
    par.route.incl = interp1(cumDistance, incl, s_x);
    par.route.max_v = interp1(cumDistance, maxSpeed, s_x);
    par.route.max_v(par.route.max_v < 51/3.6) = 60/3.6;         %initial values of max v are < 60km/h, which is the minimum velocity possible
    par.route.cumDistance = s_x;

end