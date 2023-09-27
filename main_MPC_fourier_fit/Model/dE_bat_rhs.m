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
    alpha = var(1);             % road inclination

    v_front = var(2);           % front wind velocity
    v_side = var(3);            % side wind velocity
    v_eff_front = v_front + v;
    v_eff_side = v_side;
    v_eff = sqrt(v_eff_side^2 + v_eff_front^2);
    theta_amb = var(4);
    
    G = 0;
    for i = 1:par.n_waves
        A = var((i - 1) * 3 + 5);
        B = var((i - 1) * 3 + 6);
        C = var((i - 1) * 3 + 7);
        G = G + A * sin(B * t + C);
    end

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