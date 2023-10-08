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
% "dv": dynamics state velocity

%
%%

function dv = dv_rhs(par, states, controls, var)
    
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
          

    v_eff_front = fW + v;
    v_eff = sqrt(sW^2 + v_eff_front^2);


    % longitudinal power losses
    P_a = 0.5* par.rho_a * par.Cd * par.Af * (v_eff_front)^2 * v;                 % aerodynamic loss
    P_g = par.m_tot * par.g * sin(alpha) * v;                                     % inclination loss
    P_r = par.m_tot * par.g * par.Cr * cos(alpha) * v;                            % rolling friction loss

    % electric motor traction power
    P_mot_mec = P_mot_mec_sigmoid(P_mot_el, par.e_mot, par.P_0, 0.01);        % smooth and continuous derivable    

    % velocity dynamics
    dv = (P_mot_mec + P_brake - P_a - P_g - P_r) / (par.m_tot * v);
end



function smooth_fun = P_mot_mec_sigmoid(P_mot_el, e_mot, P_0, k)
    % Sigmoid approximation for the transition step
    sigmoid = 1 ./ (1 + exp(-k * (P_mot_el - P_0)));
    
    % Approximation for y = P_mot_el * a - b if P_mot_el >= b
    part1 = P_mot_el .* e_mot - P_0;
    
    % Approximation for y = P_mot_el / a - b if P_mot_el < b
    part2 = P_mot_el ./ e_mot - P_0;
    
    % Smoothly interpolate between the two parts using the sigmoid
    smooth_fun = sigmoid .* part1 + (1 - sigmoid) .* part2;
end