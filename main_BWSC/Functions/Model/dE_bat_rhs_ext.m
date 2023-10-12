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