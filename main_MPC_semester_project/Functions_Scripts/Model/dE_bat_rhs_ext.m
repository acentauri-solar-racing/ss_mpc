% right hand side state of charge

% INPUT
% par: parameters struct
% states [n_states x 1]: states vector SX.sym('')
% controls [n_controls x 1]: controls vector SX.sym('')
% var [n_vars x 1]: variables vector SX.sym('')

function dE_bat = dE_bat_rhs_ext(par, states, controls, var)
    
    % states x and control u
    v = states(1);              % velocity
    E_bat = states(2);            % state of charge
    P_mot_el = controls(1);     % electric motor power
    
    % space dependent parameters
    alpha = var(1);             % road inclination
    G = var(2);                 % irradiation
    v_front = var(3);           % front wind velocity
    v_side = var(4);            % side wind velocity
    theta_amb = var(5);
    v_eff_front = v_front + v;
    v_eff_side = v_side;
    v_eff = sqrt(v_eff_side^2 + v_eff_front^2);

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