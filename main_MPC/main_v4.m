clear all
close all
clc

%% Define model and functions
import casadi.*
addpath("Functions_Scripts\");
addpath("Model\");
addpath("OfflineData\");
addpath("OnlineData\");

% (discretization step, horizon length, simulated distance, slack weight)
par = load_parameters(100, 10, 1000, 1e-3);

%% Define variables for optimizer CasaDi
[par, f, obj, X, U, P, S] = initialize_MPC(par);

%% Initialize constraints
[par, g_nlp, args] = initialize_constraints(par, f, X, U, P, S);

%% Create solver
solver = initialize_solver(par, obj, g_nlp, X, U, P, S);

%% Initialize Simulation
[par, OptRes] = run_simulation(par, args, f, solver, par.s_0, par.t_0);

%% Diagnostics
diagnostics = diagnostic(par.s_tot/par.s_step, OptRes.xdist);
feasible = diagnostics.alwaysFeasible

%% Print final time, plot data
final_time_s = OptRes.xx(3,end)      % [s] total time in s 
average_velocity = mean(OptRes.xx(1,:))
visualize_plot

