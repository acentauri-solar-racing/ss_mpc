% 0. setup
% 1. load parameters
% 2. initialize_mpc/nlp
% 3. initialize_constraints
% 4. initialize_solver
% 2. run_simulation
setup
%% Set main parameters

s_0 = 710000;        % initial position
s_f = 760000;       % final position
step = 100;         % simulation step

t_0 = 0;            % initial time
N = 10;             % prediction horizon

N_NLP = round((s_f-s_0)/step);

wS1 = 1e-3;         % slack 1 variable weight
wS2 = 10;           % slack 2 variable weight

% 1 run, ~= 1 not run
run_mpc = 1;
run_nlp = 1;

% initialize structs
par = get_car_param();
weather = struct;
OptRes = struct;
OptResNLP = struct;

if run_mpc == 1
    %% Load parameters
    % (par struct, s_0, t_0, discretization step, horizon length, simulated distance, slack weight)
    par = get_mpc_param(par, s_0, t_0, step, N, s_f-s_0, wS1, wS2);
    % add route parameters
    par = load_route(par);
    % add DP solution
    par = load_DP(par, 'OfflineData\Full_Race_20230607_9h_30min.mat');
    % add weather data
    % weather = load_weather(weather, "C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_semester_project\OnlineData\20230926_090342_SF\preprocess");
    weather = load_weather_benchmark(weather);


    %% Define variables for optimizer CasaDi
    [par, f, obj, X, U, P, S1, S2] = initialize_mpc(par);

    %% Initialize constraints
    [par, g_mpc, args] = initialize_constraints_mpc(par, f, X, U, P, S1, S2);
    
    %% Create solver
    solver = initialize_solver_mpc(par, obj, g_mpc, X, U, P, S1, S2);
    
    %% Initialize Simulation
    [par, OptRes] = run_simulation_mpc(par, weather, args, f, solver, par.s_0, par.t_0);
    
    %% Diagnostics
    diagnostics = diagnostic(par.s_tot/par.s_step, OptRes.xdist);
    feasible = diagnostics.alwaysFeasible
end

if run_nlp == 1
    %% NLP benchmark
    par = get_nlp_param(par, s_0, t_0, step, N_NLP, s_f-s_0, wS1, wS2);
    par = load_route(par);
    par = load_DP(par, 'OfflineData\Full_Race_20230607_9h_30min.mat');
    % weather = load_weather(weather, "C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_semester_project\OnlineData\20230926_090342_SF\preprocess");
    weather = load_weather_benchmark(weather);

    [par, f, obj, X, U, P, S1, S2] = initialize_nlp(par);
    [par, g_nlp, args] = initialize_constraints_nlp(par, f, X, U, P, S1, S2);
    solver = initialize_solver_nlp(par, obj, g_nlp, X, U, P, S1, S2);
    [par, OptResNLP] = run_simulation_nlp(par, weather, args, solver, par.s_0, par.t_0, -1, -1);
end

%% Print final time, plot data

if run_mpc == 1
    final_time_min = OptRes.xx(3,end)/60              % [min]      args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)*1.1;                     
    computational_time_MPC = OptRes.average_mpc_time*(s_f/step)
    vf_MPC = OptRes.xx(1,end)
    SoC_MPC = OptRes.xx(2,end)/par.E_bat_max
end

if run_nlp == 1
    final_time_NLP_min = OptResNLP.xx1(end,3)/60      % [min]  
    computational_time_NLP = OptResNLP.nlp_time
    vf_NLP = OptResNLP.xx1(end,1)
    SoC_NLP = OptResNLP.xx1(end,2)/par.E_bat_max
end

visualize_plot(par, OptRes, OptResNLP, run_mpc, run_nlp)
