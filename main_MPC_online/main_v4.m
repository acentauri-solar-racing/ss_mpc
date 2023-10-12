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

 addpath('..\functions\Plot');
 addpath('..\functions\MPC');
 addpath('..\functions\Model');
 addpath('..\functions\load_function');


addpath("OfflineData\");
addpath("OnlineData\");


%% initialize structs
par = get_car_param();
weather = struct;
OptRes = struct;

%% Load parameters
% (par struct, s_0, t_0, discretization step, horizon length, simulated distance, slack weight)
%par = get_mpc_param(par, init.s_0, init.t_0, init.s_step, init.N, init.N_t, init.s_f-init.s_0, init.wS1, init.wS2);
% add route parameters
par = load_route(par);
% add DP solution
par = load_DP(par, par.s_step_DP, 'OfflineData\Full_Race_20230607_9h_30min.mat');
% add weather data 
weather = load_weather(weather, par.t_0);

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
save_file(OptRes, par, weather, par.filename)