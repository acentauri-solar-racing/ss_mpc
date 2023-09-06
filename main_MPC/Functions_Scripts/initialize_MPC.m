%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 06.09.2023
%% Description 
% INPUT: 
% "par", parameters struct (see function "load parameters")

% OUTPUT : 
% "par", new parameters added,
% "f", symbolic initialization function of system dynamics to be used in
% optimization,
% "obj", objective function,
% "X", symbolic initialization of states (velocity, battery energy ,time),
% "U", symbolic initialization of control input (Electric Motor, Brake),
% "P", symbolic initialization of parameters (v0, E_bat0, t0, road inclination,
% irradiation, front wind velocity, side wind velocity, ambient
% temperature, E_bat_target),
% "S", symbolic initialization for slack variable (used for loosening SoC
% target constraint),
%%
function [par, f, obj, X, U, P, S] = initialize_MPC(par)
    %% Symbolic setup of the problem
    % Import casADi package
    import casadi.*
    
    % define symbolic variable for states
    state.v = SX.sym('state.v');
    state.E_bat = SX.sym('state.E_bat');
    state.time = SX.sym('state.time');
    
    states = [state.v; state.E_bat; state.time];
    par.n_states = length(states);
    
    % define symbolic variable for control input
    control.u = SX.sym('u');
    control.u_brake = SX.sym('u_brake');
    controls = [control.u; control.u_brake]; 
    par.n_controls = length(controls);
    
    % define symbolic variable for various parameters  
    var.alpha = SX.sym('var.alpha');                
    var.G = SX.sym('var.G');                        
    var.v_front = SX.sym('var.v_front');            
    var.v_side = SX.sym('var.v_side');             
    var.theta_amb = SX.sym('var.theta_amb');        
    variables = [var.alpha; var.G; var.v_front; var.v_side; var.theta_amb];
    par.n_vars = length(variables);
    
    % setup symbolic state dynamics right hand side 
    % define symbolic function using rhs
    % dv_rhs == longitudinal velocity
    % dE_bat_rhs == battery energy
    % 1 == time (clock)
    rhs = [dv_rhs(par,states,controls,variables); dE_bat_rhs(par,states,controls,variables); 1]; % [dv, dE_bat, dt]
    f = Function('f', {states, controls, variables},{rhs});          
    
    % setup symbolic definition of states over the entire prediction horizon
    X = SX.sym('X', par.n_states, (par.N+1));  

    % setup symbolic definition of control inputs over the entire horizon
    U = SX.sym('U', par.n_controls, par.N);     

    % setup symbolic definition of parameters over the entire horizon
    P = SX.sym('P', par.n_states + par.N*par.n_vars +1);               
    
    % setup symbolic definition of slack variable, needed only for the last
    % prediction horizon point
    S = SX.sym('S', 1, 1); 
    
    %% OBJECTIVE FUNCTION
    
    % Objective function initialization
    obj = 0; 

    % Sum of every horizon point
    for k = 1:par.N     
        % states at moment "k"
        st = X(:,k);  
        % control at moment "k"
        con = U(:,k);  
        
        % minimize time
        % 1/st(1) == 1/v == t
        obj = obj + 1/st(1);                                                 
        
    end

    % Add slack variable to objective function
    % S != 0, the SoC target constraint is loosened 
    obj = obj + par.slack_weight * S;

end