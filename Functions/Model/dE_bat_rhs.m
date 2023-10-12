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
    
    % solar irradiation
    G_1 = var(2);
    G_2 = var(3);
    G_3 = var(4);
    G = G_1*(t/60/15)^2 + G_2*(t/60/15) + G_3;

    fW_1 = var(5);
    fW_2 = var(6);
    fW_3 = var(7);
    fW = fW_1*(t/60/15)^2 + fW_2*(t/60/15) + fW_3;          

    sW_1 = var(8);
    sW_2 = var(9);
    sW_3 = var(10);
    sW = sW_1*(t/60/15)^2 + sW_2*(t/60/15) + sW_3;       

    temp_1 = var(11);
    temp_2 = var(12);
    temp_3 = var(13);
    temp = temp_1*(t/60/15)^2 + temp_2*(t/60/15) + temp_3;       
    
    %PV Model
    % eta_loss = par.eta_wire * par.eta_MPPT * par.eta_mismatch;
    
    theta_PV = temp + (par.theta_NOCT - 20) * par.S / 80;
    eta_CF = 1 - par.lambda_PV*(theta_PV - par.theta_STC);

    P_PV = par.A_PV * G * par.eta_PV_tot * eta_CF;
    
    %E_bat dynamics
    dE_bat = (P_PV - P_mot_el);          
end