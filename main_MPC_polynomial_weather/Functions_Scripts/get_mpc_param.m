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
        par.N = N;                      % [-] Horizon length
        par.N_t = par.N*par.s_step/15   % [s] Horizon length in seconds
        %par.N_t = 60*30;

        par.s_step_DP = 10000;          % [m] do not change
        
        par.s_0 = s_0;               % initial position of the simulation   
        par.s_tot =  s_tot;             % [m] simulated distance from initial position
        
        par.t_0 = t_0;              % [s], t_0 = 0 == 8.00, 0.5 == 8.30, 1 == 9.00

        %% Slack variable
        par.slack_weight = slack_weight;
        par.slack_max = inf;
        
        %% Do not change
        par.s_final = 3000000;           % [m] total distance (for parameters)
        
        %% Route and Weather 
        % Interpolated Data 
        % To be changed with updated data

        % Loading Route Data
        par.Route = load_route(par.s_step,par.s_final,par);
        par.Route.max_v(par.Route.max_v < 51/3.6) = 61/3.6;         %initial values of max v are < 60km/h, which is the minimum velocity possible
        % Loading Weather Data
        %%
        [par.timeline, par.cumdist, par.G_data, par.fW_data, par.sW_data, par.temp_data] = load_weather();

        %% DP Data 
        % Load battery target from DP
        % Should be adapted with new data
        % 300, as it was used a step resolution of 10000 m for the entire
        % race 3000000 m
        par.E_bat_target_DP_raw = load('OfflineData\Full_Race_20230607_9h_30min.mat').OptRes.states.E_bat*3600; % [Wh = W * 3600s = J]
        par.E_bat_target_DP = interp1(linspace(0,300,300+1), par.E_bat_target_DP_raw, linspace(0,300,(300)*round(par.s_step_DP/par.s_step)+1));
        
        par.v_DP_raw = load('OfflineData\Full_Race_20230607_9h_30min.mat').OptRes.states.V; % [m/s]
        par.v_DP = interp1(linspace(0,300,300+1), par.v_DP_raw, linspace(0,300,(300)*round(par.s_step_DP/par.s_step)+1));
        
        par.P_mot_el_DP_raw = load('OfflineData\Full_Race_20230607_9h_30min.mat').OptRes.inputs.P_mot_el; % [m/s]
        par.P_mot_el_DP = interp1(linspace(1,300,300), par.P_mot_el_DP_raw , linspace(1,300,300*round(par.s_step_DP/par.s_step)));
        

end