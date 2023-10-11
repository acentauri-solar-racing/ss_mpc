%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 08.09.2023
%% Description 
% This function add parameters to the struct "par" which are used in the
% MPC simulation

% INPUT: 
% "par" (struct): parameters struct
% "s_0" (int): initial position
% "t_0" (int): initial time clock
% "s_step" (int): MPC discretization step
% "N" (int): horizon length
% "s_tot" (int): total simulated distance
% "S1_weight" (int): slack variable weight in the cost objective function (SoC
% target)
% "S2_weight" (int): slack variable weight in the cost objective function (max
% velocity)

% OUTPUT : 
% "par": added mpc parameters

%% General
function par = get_mpc_param(par, s_0, t_0, s_step, N, N_t, s_tot, S1_weight, S2_weight)
        %% Discretization variables
        par.s_step = s_step;                % [m]
        par.N = N;                      % [-] Horizon length
        par.N_t = N_t;   % [s] Horizon length in seconds

        par.s_step_DP = 10000;          % [m] do not change
        
        par.s_0 = s_0;               % initial position of the simulation   
        par.s_tot =  s_tot;             % [m] simulated distance from initial position to final position
        
        par.t_0 = t_0;              % [s], t_0 = 0 == 8.00, 0.5 == 8.30, 1 == 9.00

        %% Slack variable
        par.S1_weight = S1_weight;
        par.S2_weight = S2_weight;

        par.slack_max = inf;
        
        %% Do not change
        par.s_final = 3000000;           % [m] total distance (for parameters)
        

end