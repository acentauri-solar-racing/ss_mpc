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
% "S1_weight": slack variable weight in the cost objective function (SoC
% target)
% "S2_weight": slack variable weight in the cost objective function (max
% velocity)

% OUTPUT : 
% "par": added mpc parameters

%% General
function par = get_mpc_param(par, s_0, t_0, s_step, N, s_tot, S1_weight, S2_weight)
        %% Discretization variables
        par.s_step = s_step;                % [m]
        par.N = N;                      % [-] Horizon length
        par.N_t = par.N*par.s_step/15;   % [s] Horizon length in seconds
        %par.N_t = 60*30;

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