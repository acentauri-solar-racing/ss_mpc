% 0. setup
% 1. load parameters
% 2. initialize_mpc/nlp
% 3. initialize_constraints
% 4. initialize_solver
% 2. run_simulation
clc
clear
clearvars
close all

addpath(genpath('.\..\Functions\Casadi\..'));
import casadi.*

warning('off');
addpath('..\..\ss_offline_data\route\BWSC');
addpath('..\..\ss_offline_data\parameters');
addpath('..\..\ss_offline_data\route');

addpath('..\functions\');
addpath('..\functions\Plot');
addpath('..\functions\MPC');
addpath('..\functions\Model');
addpath('..\functions\load_function');
addpath(genpath('..\functions\Casadi'));


addpath("OfflineData\");
addpath("OnlineData\");


%% initialize structs
par = get_car_param();
weather = struct;
OptRes = struct;

%% Load parameters
% add route parameters
par = load_route(par);
% add DP solution
par = load_DP(par, par.s_step_DP, 'OfflineData\Full_Race_20230607_9h_30min.mat');
% add weather data 
weather = load_weather(weather, par.t_0, 1);

%% Define variables for optimizer CasaDi
[par, mpc.f, mpc.obj, mpc.X, mpc.U, mpc.P, mpc.S1, mpc.S2] = initialize_mpc(par);

%% Initialize constraints
[par, mpc.g, mpc.args] = initialize_constraints_mpc(par, mpc.f, mpc.X, mpc.U, mpc.P, mpc.S1, mpc.S2);

%% Create solver
mpc.solver = initialize_solver_mpc(par, mpc.obj, mpc.g, mpc.X, mpc.U, mpc.P, mpc.S1, mpc.S2);

%% Initialize Simulation
[par, weather, OptRes] = run_simulation_BWSC(par, weather, mpc.args, mpc.solver, -1, -1);

%% Print final time, plot data
visualize_plot_BWSC(par, weather, OptRes)

%% Save file
save_file_BWSC(OptRes, par, weather)