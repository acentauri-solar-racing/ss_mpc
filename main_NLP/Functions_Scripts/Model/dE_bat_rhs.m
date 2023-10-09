%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% This function outputs the right hand side of the Energy Battery dynamic
% differential
% INPUT: 
% "par": parameters struct (see function "load parameters")
% "states": (velocity, battery energy, time)
% "controls": (electric motor, brake)
% "var": (inclination prediction, polynomials parameters weather fit)

% OUTPUT : 
% "dE_bat": dynamics battery energy

%% right hand side state of charge
function dE_bat = dE_bat_rhs(par, states, controls, var)
    
    % states x and control u
    v = states(1);              % velocity
    E_bat = states(2);            % state of charge
    t = states(3);
    P_mot_el = controls(1);     % electric motor power
    P_brake = controls(2);

    % road inclination
    alpha = var(1);             
    
    % solar irradiation
    G_1 = var(2);
    G_2 = var(3);
    G_3 = var(4);
    G_4 = var(5);
    
    % front wind velocity
    fW_1 = var(6);
    fW_2 = var(7);
    fW_3 = var(8);
    fW_4 = var(9);

    % air density
    rho_1 = var(10);
    rho_2 = var(11);
    rho_3 = var(12);
    rho_4 = var(13);

    % temperature
    temp_1 = var(14);
    temp_2 = var(15);
    temp_3 = var(16);
    temp_4 = var(17);

    % polynomial fit of weather variables
    G = G_1*(t/60/15)^3 + G_2*(t/60/15)^2 + G_3*(t/60/15) + G_4;
    fW = fW_1*(t/60/15)^3 + fW_2*(t/60/15)^2 + fW_3*(t/60/15) + fW_4;
    rho = rho_1*(t/60/15)^3 + rho_2*(t/60/15)^2 + rho_3*(t/60/15) + rho_4;
    temp = temp_1*(t/60/15)^3 + temp_2*(t/60/15)^2 + temp_3*(t/60/15) + temp_4;

    temp_NOCT = 50;            % [°C]
    S = 80;                     % [mW/cm^2]
    
    %PV Model
    theta_PV = temp + (temp_NOCT - 20) * S / 80;
    eta_CF = 1 - par.lambda_PV*(theta_PV - par.theta_STC);

    P_PV = par.A_PV * G * par.eta_PV_tot * eta_CF;
    
    %E_bat dynamics
    dE_bat = (P_PV - P_mot_el);          
end