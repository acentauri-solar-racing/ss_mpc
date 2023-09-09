%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 08.09.2023
%% Description 
% This function add parameters to the struct "par" which are used in the
% MPC simulation

% INPUT: 
% "par": parameters struct
% "s_0": initial position
% "t_0": initial time clock
% "s_step": MPC discretization step
% "N": horizon length
% "s_tot": total simulated distance
% "slack weight": slack variable weight in the cost objective function

% OUTPUT : 
% "par": added mpc parameters

%% General
function par = get_mpc_param(par, s_0, t_0, s_step, N, s_tot, slack_weight)
        %% Discretization variables
        par.s_step = s_step;                % [m]
        par.N = N;             % [-] Horizon length

        par.s_step_DP = 10000;          % [m] do not change
        
        par.s_0 = s_0;               % initial position of the simulation   
        par.s_tot =  s_tot;             % [m] simulated distance from initial position
        
        par.t_0 = t_0;              % [s], t_0 = 0 == 8.00, 0.5 == 8.30, 1 == 9.00

        %% Slack variable
        par.slack_weight = slack_weight;
        par.slack_max = inf;
        
        %% Do not change
        par.s_final = 3000000;           % [m] total distance (for parameters)
        par.hour_start = 8;
        par.hour_end = 17;
        
        %% Route and Weather 
        % Interpolated Data 
        % To be changed with updated data

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
        

end