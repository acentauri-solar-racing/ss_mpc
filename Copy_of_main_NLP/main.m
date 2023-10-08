% 0. setup
% 1. load parameters
% 2. initialize_mpc/nlp
% 3. initialize_constraints
% 4. initialize_solver
% 2. run_simulation
setup
%% Set main parameters
v_0 = -1;
SoC_0 = -1;
s_0 = 0;        % initial position
t_0 = 0;            % initial time

s_f = 10000;       % final position

step = 100;         % simulation step
N_NLP = round((s_f-s_0)/step);


wS1 = 1e-3;         % slack 1 variable weight
wS2 = 1e4;           % slack 2 variable weight
%%
% initialize structs
par = get_car_param();
weather = struct;
OptResNLP = struct;

%% NLP benchmark
par = get_nlp_param(par, s_0, t_0, step, N_NLP, s_f-s_0, wS1, wS2);
par = load_route(par);
par = load_DP(par, 'OfflineData\Full_Race_20230607_9h_30min.mat');
weather = load_weather(weather, "C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_semester_project\OnlineData\20230926_090342_SF\preprocess\");
%weather = load_weather_benchmark(weather);

[par, f, obj, X, U, P, S1, S2] = initialize_nlp(par);
[par, g_nlp, args] = initialize_constraints_nlp(par, f, X, U, P, S1, S2);
solver = initialize_solver_nlp(par, obj, g_nlp, X, U, P, S1, S2);
[par, OptResNLP] = run_simulation_nlp(par, weather, args, solver, par.s_0, par.t_0, -1, 0.99, -1, -1);

%% Print final time, plot data

final_time_NLP_min = OptResNLP.xx1(end,3)/60      % [min]  
computational_time_NLP = OptResNLP.nlp_time
vf_NLP = OptResNLP.xx1(end,1)
SoC_NLP = OptResNLP.xx1(end,2)/par.E_bat_max

visualize_plot(par, OptResNLP)
