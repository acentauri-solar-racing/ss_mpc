% 0. setup
% 1. load parameters
% 2. initialize_mpc/nlp
% 3. initialize_constraints
% 4. initialize_solver
% 2. run_simulation
setup


%% Set main parameters
t_0 = 60*60*8.5;                          % initial time in s
%t_0 = get_machine_time_s(); % machine time
s_0 = 0;                    % [m] initial position, s_0 >= DP_s_0

s_f = 10000;        % [m] final position, s_f > s_0
step = 100;          % [m] simulation step
N = 100;             % [-] prediction horizon
N_t = 60*15*2.5;        % [-] weather fit horizon

format = 'yyyymmdd_HHMMSS';
filename = [datestr(datetime,format)+"_"+num2str(s_0)+"_"+num2str(s_f)];
%% DP settings
DP_s_0 = 0;    % [m] DP initial position
DP_step = 10000;    % [m] DP resolution
%% Slack Variable weights
wS1 = 1e-3;          % slack 1 variable weight
wS2 = 1e4;           % slack 2 variable weight

%% initialize structs
par = get_car_param();
weather = struct;
OptRes = struct;
OptResNLP = struct;

%% Load parameters
% (par struct, s_0, t_0, discretization step, horizon length, simulated distance, slack weight)
par = get_mpc_param(par, s_0, t_0, step, N, N_t, s_f-s_0, wS1, wS2);
% add route parameters
par = load_route(par);
% add DP solution
par = load_DP(par, DP_step, 'OfflineData\Full_Race_20230607_9h_30min.mat');
% add weather data 
weather = load_weather(weather, par.t_0, "G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\Forecast\20231008_094711_SF\preprocess\");
%weather = load_weather_benchmark(weather);


%% Define variables for optimizer CasaDi
[par, f, obj, X, U, P, S1, S2] = initialize_mpc(par);

%% Initialize constraints
[par, g_mpc, args] = initialize_constraints_mpc(par, f, X, U, P, S1, S2);

%% Create solver
solver = initialize_solver_mpc(par, obj, g_mpc, X, U, P, S1, S2);

%% Initialize Simulation
[par, OptRes] = run_simulation_mpc_online(par, weather, args, f, solver, par.s_0, DP_s_0, par.t_0);

%% Diagnostics
% diagnostics = diagnostic(par.s_tot/par.s_step, OptRes.xdist);
% feasible = diagnostics.alwaysFeasible

%% Print final time, plot data

final_time_min = (OptRes.xx(3,end)-t_0)/60              % [min]      args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)*1.1;                     
computational_time_MPC = OptRes.average_mpc_time/((s_f-s_0)/step)
OptRes.average_mpc_time = OptRes.average_mpc_time
vf_MPC = OptRes.xx(1,end)
SoC_MPC = OptRes.xx(2,end)/par.E_bat_max


visualize_plot(par, OptRes)
%% Save file
save('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\'+filename, "OptRes","par","weather")
writematrix(OptRes.xx,'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\'+filename+'_state.csv')