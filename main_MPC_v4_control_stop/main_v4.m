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
addpath("Functions_Scripts\load_function");
addpath("Functions_Scripts\MPC");
addpath("Functions_Scripts\NLP");


addpath("Model\");
addpath("OfflineData\");
addpath("OnlineData\");
%%
s_0 = 0;
t_0 = 0;
step = 100; 
s_f = 330000;
N = 20; 
N_NLP = round(s_f/step);
wS1 = 1e-3;
wS2 = 10;

%% Load parameters

% common car parameters
par = get_car_param();
% (par struct, s_0, t_0, discretization step, horizon length, simulated distance, slack weight)
par = get_mpc_param(par, s_0, t_0, step, N, s_f, wS1, wS2);

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
feasible = diagnostics.alwaysFeasible;


%% NLP benchmark
% par.N_NLP = round(par.s_tot/par.s_step);
% 
% par = get_mpc_param(par, s_0, t_0, step, N_NLP, s_f, wS1, wS2);
% 
% [par, f, obj, X, U, P, S1, S2] = initialize_MPC(par);
% 
% [par, g_nlp, args] = initialize_constraints_NLP(par, f, X, U, P, S1, S2);
% 
% solver = initialize_solver_NLP(par, obj, g_nlp, X, U, P, S1, S2);
% 
% [par, OptResNLP] = run_simulation_nlp(par, args, solver, par.s_0, par.t_0, par.final_velocity, par.final_E_bat);
% %% Print final time, plot data
% 
% final_time_min = OptRes.xx(3,end)/60              % [min]      args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)*1.1;                     
% final_time_NLP_min = OptResNLP.xx1(end,3)/60      % [min]  
% 
% computational_time_NLP = OptResNLP.nlp_time
% computational_time_MPC = OptRes.average_mpc_time*(s_f/step)
% 
% vf_MPC = OptRes.xx(1,end)
% vf_NLP = OptResNLP.xx1(end,1)
% 
% SoC_MPC = OptRes.xx(2,end)/par.E_bat_max
% SoC_NLP = OptResNLP.xx1(end,2)/par.E_bat_max
% 
visualize_plot
