%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% INPUT: 
% "par": parameters struct (see function "load parameters")
% "states": (velocity, battery energy, time)
% "controls": (electric motor, brake)
% "var": (inclination prediction, polynomials parameters weather fit)

% OUTPUT : 
% "dE_bat": dynamics battery energy
%%

function dE_bat = dE_bat_rhs_ext(par, states, controls, var)
    
    % states x and control u
    v = states(1);              % velocity
    E_bat = states(2);            % state of charge
    t = states(3);
    P_mot_el = controls(1);     % electric motor power
    P_brake = controls(2);

    % space dependent parameters

    % road inclination
    alpha = var(1);             
    
    % solar irradiation
    G_1 = var(2);
    G_2 = var(3);
    G_3 = var(4);
    G_4 = var(5);

    G = G_1*(t/60/15)^3 + G_2*(t/60/15)^2 + G_3*(t/60/15) + G_4;

    fW_1 = var(6);
    fW_2 = var(7);
    fW_3 = var(8);
    fW_4 = var(9);

    fW = fW_1*(t/60/15)^3 + fW_2*(t/60/15)^2 + fW_3*(t/60/15) + fW_4;

    sW_1 = var(10);
    sW_2 = var(11);
    sW_3 = var(12);
    sW_4 = var(13);

    sW = sW_1*(t/60/15)^3 + sW_2*(t/60/15)^2 + sW_3*(t/60/15) + sW_4;

    temp_1 = var(14);
    temp_2 = var(15);
    temp_3 = var(16);
    temp_4 = var(17);

    temp = temp_1*(t/60/15)^3 + temp_2*(t/60/15)^2 + temp_3*(t/60/15) + temp_4;
    
    theta_NOCT = 50;             % [°C]
    S = 80;                     % [mW/cm^2]

    %PV Model
    eta_loss = par.eta_wire * par.eta_MPPT * par.eta_mismatch;
    
    theta_PV = theta_amb + (theta_NOCT - 20) * S / 80; % [°C]
    eta_CF = 1 - par.lambda_PV*(theta_PV - par.theta_STC);

    P_PV = par.A_PV * G * par.eta_PV * eta_CF * eta_loss;
    
    % Battery Power
    P_bat = (P_mot_el - P_PV);     
    
    % Battery Current
    I_bat = (par.U_bat_oc - sqrt(par.U_bat_oc^2 - 4 * par.R_bat * P_bat))/(2*par.R_bat);
    
    % Charge flow
    Q_bat_dot = -par.eta_coul * I_bat;

    % dE_bat
    dE_bat = Q_bat_dot * par.U_bat_oc;

end