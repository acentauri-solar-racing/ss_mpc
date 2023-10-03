clear all
close all
clc

%% Define model and functions
import casadi.*

warning('off');
addpath('..\..\ss_offline_data\route\BWSC');

addpath('..\..\ss_offline_data\parameters');
addpath('..\..\ss_offline_data\route');

addpath("Functions_Scripts\");
addpath("Functions_Scripts\load_function\");

addpath("Model\");
addpath("OfflineData\");
addpath("OnlineData\");

%% Load parameters

% common car parameters
par = get_car_param();
% (par struct, s_0, t_0, discretization step, horizon length, simulated distance, slack weight)
par = get_mpc_param(par, 450000, 5000, 100, 20, 50000, 1e-3, 10);

%% Define variables for optimizer CasaDi
[par, f, obj, X, U, P, S1, S2] = initialize_MPC(par);

%% Initialize constraints
[par, g_nlp, args] = initialize_constraints(par, f, X, U, P, S1, S2);

%% Create solver
solver = initialize_solver(par, obj, g_nlp, X, U, P, S1, S2);

%% Initialize Simulation
[par, OptRes] = run_simulation_mpc(par, args, f, solver, par.s_0, par.t_0);

%% Diagnostics
diagnostics = diagnostic(par.s_tot/par.s_step, OptRes.xdist);
feasible = diagnostics.alwaysFeasible

%% Print final time, plot data
final_time_s = OptRes.xx(3,end)      % [s] total time in s         args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)*1.1;                     

average_velocity = mean(OptRes.xx(1,:))
visualize_plot
