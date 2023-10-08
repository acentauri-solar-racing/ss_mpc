%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% This function initializes the symbolic variable and functions using
% casADi that are used in the solver optimization and throughout the entire
% code.
% To learn more about MPC and casADi it's HIGHLY suggested to watch this
% tutorial video 
% https://www.youtube.com/watch?v=RrnkPrcpyEA

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
% "S1"/"S2", symbolic initialization for slack variable (used for loosening SoC
% target constraint),
%%
function [par, f, obj, X, U, P, S1, S2] = initialize_nlp(par)
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

    var.G_1 = SX.sym('var.G_1');
    var.G_2 = SX.sym('var.G_2');                        
    var.G_3 = SX.sym('var.G_3');
    var.G_4 = SX.sym('var.G_4');                        

    var.fW_1 = SX.sym('var.fW_1');
    var.fW_2 = SX.sym('var.fW_2');                        
    var.fW_3 = SX.sym('var.fW_3'); 
    var.fW_4 = SX.sym('var.fW_4'); 

    var.sW_1 = SX.sym('var.sW_1');
    var.sW_2 = SX.sym('var.sW_2');                        
    var.sW_3 = SX.sym('var.sW_3');    
    var.sW_4 = SX.sym('var.sW_4');    

    var.temp_1 = SX.sym('var.temp_1');
    var.temp_2 = SX.sym('var.temp_2');                        
    var.temp_3 = SX.sym('var.temp_3');    
    var.temp_4 = SX.sym('var.temp_4');    

    variables = [var.alpha; 
                var.G_1; var.G_2; var.G_3; var.G_4;
                var.fW_1; var.fW_2; var.fW_3; var.fW_4;
                var.sW_1; var.sW_2; var.sW_3;  var.sW_4;
                var.temp_1; var.temp_2; var.temp_3; var.temp_4];

    par.n_vars = length(variables);
    
    % setup symbolic state dynamics right hand side 
    % define symbolic function using rhs
    % dv_rhs == longitudinal velocity
    % dE_bat_rhs == battery energy
    % 1 == time (clock)
    rhs = [dv_rhs(par,states,controls,variables); dE_bat_rhs_ext(par,states,controls,variables); 1]; % [dv, dE_bat, dt]
    f = Function('f', {states, controls, variables},{rhs});          
    
    % setup symbolic definition of states over the entire prediction horizon
    X = SX.sym('X', par.n_states, (par.N+1));  

    % setup symbolic definition of control inputs over the entire horizon
    U = SX.sym('U', par.n_controls, par.N);     

    % setup symbolic definition of parameters over the entire horizon
    P = SX.sym('P', par.n_states + par.N + 16 +1);               
    
    % setup symbolic definition of slack variable, needed only for the last
    % prediction horizon point
    S1 = SX.sym('S1', 1, 1); 
    S2 = SX.sym('S2', 1, 1); 

    
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
    obj = obj + par.S1_weight * S1 + par.S2_weight * S2;

end