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
addpath("Functions_Scripts\Model");
addpath("Functions_Scripts\load_function");
addpath("Functions_Scripts\MPC");
addpath("Functions_Scripts\NLP");

addpath("OfflineData\");
addpath("OnlineData\");
%% Set main parameters

s_0 = 310000;        % initial position
s_f = 330000;       % final position
step = 100;         % simulation step

t_0 = 60*60*3;            % initial time
N = 10;             % prediction horizon

N_NLP = round((s_f-s_0)/step);

wS1 = 1e-3;         % slack 1 variable weight
wS2 = 10;           % slack 2 variable weight

run_mpc = 1;
run_nlp = 1;

if run_mpc == 1;
    %% Load parameters
    
    % common car parameters
    par = get_car_param();
    % (par struct, s_0, t_0, discretization step, horizon length, simulated distance, slack weight)
    par = get_mpc_param(par, s_0, t_0, step, N, s_f-s_0, wS1, wS2);
    
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
end

if run_nlp == 1
    %% NLP benchmark
    par = get_car_param();
    par = get_mpc_param(par, s_0, t_0, step, N_NLP, s_f-s_0, wS1, wS2);
    
    [par, f, obj, X, U, P, S1, S2] = initialize_MPC(par);
    
    [par, g_nlp, args] = initialize_constraints_NLP(par, f, X, U, P, S1, S2);
    
    solver = initialize_solver_NLP(par, obj, g_nlp, X, U, P, S1, S2);
    
    [par, OptResNLP] = run_simulation_nlp(par, args, solver, par.s_0, par.t_0, -1, -1);
end

%% Print final time, plot data

if run_mpc == 1;
    final_time_min = OptRes.xx(3,end)/60              % [min]      args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)*1.1;                     
    computational_time_MPC = OptRes.average_mpc_time*(s_f/step)
    vf_MPC = OptRes.xx(1,end)
    SoC_MPC = OptRes.xx(2,end)/par.E_bat_max
end

if run_nlp == 1;
    final_time_NLP_min = OptResNLP.xx1(end,3)/60      % [min]  
    computational_time_NLP = OptResNLP.nlp_time
    vf_NLP = OptResNLP.xx1(end,1)
    SoC_NLP = OptResNLP.xx1(end,2)/par.E_bat_max
end

visualize_plot
