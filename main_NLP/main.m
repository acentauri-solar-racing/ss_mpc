% 0. setup
% 1. load parameters
% 2. initialize_mpc/nlp
% 3. initialize_constraints
% 4. initialize_solver
% 2. run_simulation
setup
%% Set main parameters
v_0 = -1;
SoC_0 = 0.70;
v_target = -1;
SoC_target = -1;      % -1

s_0 = 500000;        % initial position
t_0 = 10000; %60*15*0;            % initial time from which the data starts

s_f = s_0 + 50000;       % final position

step = 100;         % simulation step
N_NLP = round((s_f-s_0)/step);
N_t = 60*15*2; % [s] prediction weather fit

wS1 = 1e-3;         % slack 1 variable weight
wS2 = 1e4;           % slack 2 variable weight
%%
% initialize structs
par = get_car_param();
weather = struct;
OptResNLP = struct;

%% NLP benchmark
par = get_nlp_param(par, s_0, t_0, step, N_NLP, N_t, s_f-s_0, wS1, wS2);
par = load_route(par);
par = load_DP(par, 'OfflineData\Full_Race_20230607_9h_30min.mat');
weather = load_weather(weather, s_0, "G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\Forecast\20231008_094711_SF\preprocess\");
%weather = load_weather_benchmark(weather);

[par, f, obj, X, U, P, S1, S2] = initialize_nlp(par);
[par, g_nlp, args] = initialize_constraints_nlp(par, f, X, U, P, S1, S2);
solver = initialize_solver_nlp(par, obj, g_nlp, X, U, P, S1, S2);
[par, OptResNLP] = run_simulation_nlp(par, weather, args, solver, par.s_0, par.t_0, v_0, SoC_0, v_target, SoC_target);

%% Print final time, plot data

final_time_NLP_min = (OptResNLP.xx1(end,3)-t_0)/60      % [min]  
computational_time_NLP = OptResNLP.nlp_time
vf_NLP = OptResNLP.xx1(end,1)
SoC_NLP = OptResNLP.xx1(end,2)/par.E_bat_max

visualize_plot(par, weather, OptResNLP)

% 50 km -> 15 min x 2
% 100 km -> 15 min x3/ x4