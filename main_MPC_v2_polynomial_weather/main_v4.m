clear all
close all
clc

%% Define model and functions
import casadi.*

warning('off');
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
par = get_mpc_param(par, 0, 0, 100, 10, 200000, 1e-3);

%% Define variables for optimizer CasaDi
[par, f, obj, X, U, P, S] = initialize_MPC(par);

%% Initialize constraints
[par, g_nlp, args] = initialize_constraints(par, f, X, U, P, S);

%% Create solver
solver = initialize_solver(par, obj, g_nlp, X, U, P, S);

%% Initialize Simulation
[par, OptRes] = run_simulation_mpc(par, args, f, solver, par.s_0, par.t_0);

%% Diagnostics
diagnostics = diagnostic(par.s_tot/par.s_step, OptRes.xdist);
feasible = diagnostics.alwaysFeasible

%% Print final time, plot data
final_time_s = OptRes.xx(3,end)      % [s] total time in s 
average_velocity = mean(OptRes.xx(1,:))
visualize_plot

%% NLP benchmark
% par.N_NLP = round(par.s_tot/par.s_step);
% 
% par = get_mpc_param(par, 0, 0, 100, par.N_NLP, s_tot, 1e-3);
% 
% [par, f, obj, X, U, P, S] = initialize_MPC(par);
% 
% [par, g_nlp, args] = initialize_constraints(par, f, X, U, P, S);
% 
% solver = initialize_solver(par, obj, g_nlp, X, U, P, S);
% 
% [par, OptResNLP] = run_simulation_nlp(par, args, f, solver, par.s_0, par.t_0, par.final_velocity, par.final_E_bat);
% %final_time_NLP = OptResNLP.xx(3,end)      % [s] total time in s 
