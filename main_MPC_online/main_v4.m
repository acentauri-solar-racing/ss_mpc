% 0. setup
% 1. load parameters
% 2. initialize_mpc/nlp
% 3. initialize_constraints
% 4. initialize_solver
% 2. run_simulation
clear all
close all
clc

import casadi.*

warning('off');
addpath('..\..\ss_offline_data\route\BWSC');

addpath('..\..\ss_offline_data\parameters');
addpath('..\..\ss_offline_data\route');

addpath("Functions\");
addpath("Functions\Model");
addpath("Functions\Plot");
addpath("Functions\Load_Function");
addpath("Functions\MPC");
addpath("Functions\NLP");

addpath("OfflineData\");
addpath("OnlineData\");


%% Set main parameters
% "init" struct, just for "order"; init means initialized, because here you
% initialize the data; 
init.t_0 = 60*60*15;
%init.t_0 = get_machine_time_s();   % [s] machine time (REMEMBER TO ADD TIME IF YOU ARE IN THE NIGHT)
init.s_0 = 1000000;                 % [m] initial position, s_0 >= DP_s_0

init.s_f = 1000 + init.s_0;   % [m] final position, s_f > s_0
init.s_step = 100;           % [m] simulation step
init.N = 100;                % [-] prediction horizon
init.N_t = 60*15*2.5;        % [-] weather fit horizon

init.format = 'yyyymmdd_HHMMSS';
filename.date = [datestr(datetime,init.format)+"_"+num2str(init.s_0)+"_"+num2str(init.s_f)];
%% DP settings
init.DP_s_0 = 0;    % [m] DP initial position
init.DP_step = 10000;    % [m] DP resolution
%% Slack Variable weights
init.wS1 = 1e-3;          % slack 1 variable weight
init.wS2 = 1e4;           % slack 2 variable weight

%% initialize structs
par = get_car_param();
weather = struct;
OptRes = struct;

%% Load parameters
% (par struct, s_0, t_0, discretization step, horizon length, simulated distance, slack weight)
par = get_mpc_param(par, init.s_0, init.t_0, init.s_step, init.N, init.N_t, init.s_f-init.s_0, init.wS1, init.wS2);
% add route parameters
par = load_route(par);
% add DP solution
par = load_DP(par, init.DP_step, 'OfflineData\Full_Race_20230607_9h_30min.mat');
% add weather data 
filename.folderweather = get_latest_weather('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\Forecast\');
weather = load_weather(weather, par.t_0, filename.folderweather);

%% Define variables for optimizer CasaDi
[par, mpc.f, mpc.obj, mpc.X, mpc.U, mpc.P, mpc.S1, mpc.S2] = initialize_mpc(par);

%% Initialize constraints
[par, mpc.g, mpc.args] = initialize_constraints_mpc(par, mpc.f, mpc.X, mpc.U, mpc.P, mpc.S1, mpc.S2);

%% Create solver
mpc.solver = initialize_solver_mpc(par, mpc.obj, mpc.g, mpc.X, mpc.U, mpc.P, mpc.S1, mpc.S2);

%% Initialize Simulation
[par, OptRes] = run_simulation_mpc_online(par, weather, mpc.args, mpc.f, mpc.solver, par.s_0, par.t_0);

%% Print final time, plot data
visualize_plot(par, OptRes)

%% Save file
save_file(OptRes, par, weather, filename.date)