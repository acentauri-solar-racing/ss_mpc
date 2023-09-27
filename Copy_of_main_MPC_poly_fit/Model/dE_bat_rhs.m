%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 06.09.2023
%% Description 
% INPUT: 
% "par": parameters struct (see function "load parameters")
% "states":
% "controls":
% "var":

% OUTPUT : 
% "dE_bat"

%% right hand side state of charge
function dE_bat = dE_bat_rhs(par, states, controls, var)
    
    % states x and control u
    v = states(1);              % velocity
    E_bat = states(2);            % state of charge
    t = states(3);

    P_mot_el = controls(1);     % electric motor power
    
    % space dependent parameters

    % road inclination
    alpha = var(1);             
    
    % front wind velocity
    v_front = var(2);           
    
    % side wind velocity
    v_side = var(3);            

    v_eff_front = v_front + v;
    v_eff_side = v_side;
    v_eff = sqrt(v_eff_side^2 + v_eff_front^2);

    theta_amb = var(4);
    
    % solar irradiation
    G_1 = var(5);
    G_2 = var(6);
    G_3 = var(7);
    G = G_1*(t/60)^2 + G_2*(t/60) + G_3;
    
    theta_NOCT = 50;             % [°C]
    S = 80;                     % [mW/cm^2]
    
    %PV Model
    % eta_loss = par.eta_wire * par.eta_MPPT * par.eta_mismatch;
    
    theta_PV = theta_amb + (theta_NOCT - 20) * S / 80;
    eta_CF = 1 - par.lambda_PV*(theta_PV - par.theta_STC);

    P_PV = par.A_PV * G * par.eta_PV_tot * eta_CF;
    
    %E_bat dynamics
    dE_bat = (P_PV - P_mot_el);          
end