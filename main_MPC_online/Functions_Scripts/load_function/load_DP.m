%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% This function outputs the right hand side of the velocity dynamic
% differential
% equation
% INPUT: 
% "par": parameters struct (see function "load parameters")
% "DP_step": initial position of the DP in m;
% "filenamepath": name of the path of the DP file solution; IMPORTANT THE
% DP SOLUTION CAN BE START WHEREVER, BUT HAS TO FINISH AT 3000 KM! (END OF
% THE RACE)

% OUTPUT : 
% "par": updated struct that includes DP solution

%%
function par = load_DP(par, DP_step, filenamepath)
        % end of the race
        s_tot = 3000000;
        % data length of entire race
        full_len = s_tot/par.s_step+1;
        % load DP solution
        DP_sol = load(filenamepath);
        
        % raw data from DP
        par.E_bat_target_DP_raw = load(filenamepath).OptRes.states.E_bat*3600; % [Wh = W * 3600s = J]
        par.v_DP_raw = load(filenamepath).OptRes.states.V; % [m/s]
        par.P_mot_el_DP_raw = load(filenamepath).OptRes.inputs.P_mot_el; % [W]



        % length raw data (<= data length of entire race)
        len = length(par.E_bat_target_DP_raw);
        len_u = length(par.P_mot_el_DP_raw);

        % interpolate to have the same resolution of MPC
        par.E_bat_target_DP = interp1(linspace(0,len-1,len), par.E_bat_target_DP_raw, linspace(0,len-1,(len-1)*round(DP_step/par.s_step)+1));
        par.v_DP = interp1(linspace(0,len-1,len), par.v_DP_raw, linspace(0,len-1,(len-1)*round(DP_step/par.s_step)+1)); 
        par.P_mot_el_DP = interp1(linspace(1,len_u,len_u), par.P_mot_el_DP_raw , linspace(1,len_u,len_u*round(DP_step/par.s_step)));
       
        % Create a vector filled with -1 of length (3000000 - n)
        filler = -1 * ones(1, full_len - length(par.E_bat_target_DP));
        filler_u = -1 * ones(1, full_len-1 - length(par.E_bat_target_DP));

        % Concatenate the filler vector with the original_vector, so that
        % it has the full data length of the race (full_len)
        par.E_bat_target_DP_raw = [filler, par.E_bat_target_DP_raw];
        par.v_DP_raw = [filler, par.E_bat_target_DP_raw];
        par.P_mot_el_DP_raw = [filler_u, par.P_mot_el_DP_raw];

end