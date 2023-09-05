%parameters file

%% General
par.s_step = 50;                % [m]
par.s_step_DP = 10000;          % [m] do not change

par.s_initial = 0;               % initial position of the simulation   
par.s_tot =  500000;             % [m] simulated distance from initial position
par.N_horizon = 20;             % [-] Horizon length

par.t_0 = 0 *60*60;              % [s], t_0 = 0 == 8.00, 0.5 == 8.30, 1 == 9.00

par.s_final = 3000000;           % [m] total distance (for parameters)
par.hour_start = 8;
par.hour_end = 17;

par.v_max = 120/3.6;             % abs max velocity
par.v_min = 60/3.6;

%% Longitudinal Vehicle Dynamics
par.rho_a = 1.17;           % [kg/m^3] air density 
par.Af = 0.8;               % [m^2] frontal area 
par.Cd = 0.09;              % [-] aero drag coeff 
par.Cr = 0.003;             % [-] roll fric coeff 
par.g = 9.81;               % [m/s^2] gravitational acc
par.r_w = 0.2785;           % [m] wheel radius 
par.N_f = 4;                % [-] #front bearings
par.T_f = 0.0550;           % [Nm] friction torque in one front bearing
par.N_r = 1;                % [-] #back bearings
par.T_r = 0.15;             % [Nm] friction torque back bearing
par.m_car = 150;            % [kg] car mass
par.m_driver = 80;          % [kg] driver mass
par.Theta_rot = 1.1343;     % [kgm^2] Moment of Inertia rotating parts
par.m_dt = par.Theta_rot/par.r_w^2; % [kg] mass of rotating parts
par.m_tot = par.m_car+par.m_driver+par.m_dt; % [kg] total mass

%% Electric Motor
par.gamma_gb = 1;           % [-] transmission gear box
par.e_mot = 0.97;           % [-] Motor efficiency
par.P_0 = 30;               % [W] Idle losses
par.P_el_max = 5000;        % [W] Maximal electric power
par.P_el_min = -5000;       % [W] Minimal electric power
par.T_mot_max = 45;         % [Nm] Maximal Torque
par.T_mot_min = -45;        % [Nm] Minimal Torque

%Photvoltaic Module
par.nu = -3.47;             % [-]
par.kappa = -0.0594;        % [s/m]
par.A_PV = 4;               % [-] solar panel area
par.eta_PV = 0.244;         % [-] solar panel efficiency
par.eta_wire = 0.98;        % [-] wiring efficiency
par.eta_MPPT = 0.99;        % [-] MPPT efficiency
par.eta_mismatch = 0.98;    % [-] Mismatch efficiency
par.theta_STC = 25;         % [°C] Standard Condition Temperature
par.G_0 = 1000;             % [W/m^2] Reference Global Irradiance
par.lambda_PV = 0.0029;     % [1/K] Power loss coefficient

%% Battery Pack
par.U_bat_oc = 126;         % [V] open circuit voltage
par.R_bat = 0.075;          % [Ohm] Internal Resistance
par.E_bat_max = 5200*3600;  % [W*s = J] Maximal energy capacity
par.I_bat_max = 78.4;       % [A] max current
par.I_bat_min = -39.2;      % [A] min current
par.P_bat_max = 9878.4;     % [W] max power
par.P_bat_min = -4939.4;    % [W] min power
par.SoC_max = 1;            % [-] max safe SoC
par.SoC_min = 0.1;          % [-] min safe SoC
par.eta_coul = 1;           % [-] Coulumbic efficiency

%% Route and Weather 
% Loading Route Data
par.Route = load_route(par.s_step,par.s_final,par);
par.Route.max_v(par.Route.max_v < 51/3.6) = 61/3.6;         %initial values of max v are < 60km/h, which is the minimum velocity possible
% Loading Weather Data
% Irradiance 1D time interpolation
par.G.rawdata = load('OnlineData\WeatherIrradiance.mat').irradiance.Gtotal(par.hour_start*60+1:par.hour_end*60);
par.G.t = linspace(0,60*(par.hour_end-par.hour_start),60*(par.hour_end-par.hour_start));        %vector time in MINUTES
par.G_int = griddedInterpolant(par.G.t,par.G.rawdata);                                               %remember to convert seconds into minutes!

% Temperature 2D spatial-time interpolation
par.theta.rawdata = load("OnlineData\WeatherData.mat").temperature.tempMean(:,5:14); % from 8.30 to 17.30
[par.theta.dist, par.theta.t] = ndgrid(load("OnlineData\WeatherData.mat").temperature.dist,linspace(0,9,10));
par.theta_int = griddedInterpolant(par.theta.dist,par.theta.t,par.theta.rawdata);

% Front Wind 2D spatial-time interpolation
par.wind_front.rawdata = load("OnlineData\WeatherData.mat").wind.frontWind(:,5:14); % from 8.30 to 17.30
[par.wind.dist, par.wind.t] = ndgrid(load("OnlineData\WeatherData.mat").wind.dist(1:end-1),linspace(0,9,10));
par.wind_front_int = griddedInterpolant(par.wind.dist,par.wind.t,par.wind_front.rawdata);

% Side Wind 2D spatial-time interpolation
par.wind_side.rawdata = load("OnlineData\WeatherData.mat").wind.sideWind(:,5:14); % from 8.30 to 17.30
par.wind_side_int = griddedInterpolant(par.wind.dist,par.wind.t,par.wind_side.rawdata);
%% DP Data 
% Load battery target from DP
par.E_bat_target_DP_raw = load('OfflineData\Full_Race_20230607_9h_30min.mat').OptRes.states.E_bat*3600; % [Wh = W * 3600s = J]
par.E_bat_target_DP = interp1(linspace(0,300,300+1), par.E_bat_target_DP_raw, linspace(0,300,(300+1)*round(par.s_step_DP/par.s_step)));

par.v_DP_raw = load('OfflineData\Full_Race_20230607_9h_30min.mat').OptRes.states.V; % [m/s]
par.v_DP = interp1(linspace(0,300,300+1), par.v_DP_raw, linspace(0,300,(300+1)*round(par.s_step_DP/par.s_step)));

par.P_mot_el_DP_raw = load('OfflineData\Full_Race_20230607_9h_30min.mat').OptRes.inputs.P_mot_el; % [m/s]
par.P_mot_el_DP = interp1(linspace(1,300,300), par.P_mot_el_DP_raw , linspace(1,300,300*round(par.s_step_DP/par.s_step)));

%% Slack variable

par.slack_weight = 1e-3;
par.slack_max = inf;