%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% This function initializes the constraints used in the nlp.
% It first defines the symbolic constraint functions in the vector "g_nlp"
% and set later the lower/upper boundaries in the "args" struct, of both
% constraints and state/control variables

% args.lbg <= g_nlp <= args.ubg
% args.lbx <= opt_variables <= args.ubx

% INPUT: 
% "par": parameters struct,
% "f": symbolic initialization function of system dynamics to be used in
% optimization,
% "X": symbolic initialization of states (velocity, battery energy ,time),
% "U": symbolic initialization of control input (Electric Motor, Brake),
% "P": symbolic initialization of parameters (v0, E_bat0, t0, road inclination,
% irradiation, front wind velocity, side wind velocity, ambient
% temperature, E_bat_target),
% "S1", "S2": symbolic initialization for slack variable (used for loosening SoC
% target constraint),

% OUTPUT : 
% "par": new parameters added,
% "g_nlp": constraints vector,
%          "g,x0": constraint that enforces initial conditions,
%          "g.st_next": constraints of multiple shooting, where it is
%          enforced x(k+1) = f(x(k),u(k),w(k))

% "args": struct used in optimizer, variable that can be updated during
%         simulation,
%         "args.lbg": contains lower bounds for constraints vector "g_nlp"
%         "args.ubg": contains upper bounds for constraint vector "g_nlp"
%         "args.lbx": contains lower bounds for optimization variables vector 
%         "args.ubx": contains upper bounds for optimization variables vector 
%          N.B. optimization variables vector contains states prediction,
%          control input prediction, slack variable

%%
function [par, g_nlp, args] = initialize_constraints_nlp(par, f, X, U, P, S1, S2)
    %% Define constraints vector
    % initialize constraint vector
    g_nlp = [];                                             
    
    % initial state x(k) = x(0)
    st  = X(:,1); 
    
    % initial condition constraints
    g.x0 = st-P(1:3); 
    
    % initialize multiple shooting constraint
    g.st_next = [];
    
    % constraint multiple shooting 
    for k = 1:par.N     
        % states at step k, variable to optimize x_opt(k)
        st = X(:,k);  
        % controls at step k, variable to optimize u_opt(k)
        con = U(:,k); 

        % define inclination/weather conditions at timestep k
        w = [P(par.n_states+ k);                               % inclination    
             P(par.n_states+par.N + 1);                        % G1
             P(par.n_states+par.N + 2);                        % G2
             P(par.n_states+par.N + 3);                        % G3
             P(par.n_states+par.N + 4);                        % fW1
             P(par.n_states+par.N + 5);                        % fW2
             P(par.n_states+par.N + 6);                        % fW3
             P(par.n_states+par.N + 7);                        % sW1
             P(par.n_states+par.N + 8);                        % sW2
             P(par.n_states+par.N + 9);                        % sW3
             P(par.n_states+par.N + 10);                       % temp1
             P(par.n_states+par.N + 11);                       % temp2
             P(par.n_states+par.N + 12);                       % temp3
             P(par.n_states+par.N + 13);                        % sW3
             P(par.n_states+par.N + 14);                       % temp1
             P(par.n_states+par.N + 15);                       % temp2
             P(par.n_states+par.N + 16);                       % temp3
             ];         
            
        % Runge Kutta 4th order, compute x(k+1) = f(x_opt(k),u_opt(k),w(k))
        % NB par.s_step/st(1) == ds/v == dt

        % RK weights
        k1 = f(st, con, w);                                                    
        k2 = f(st + (par.s_step/st(1))/2*k1, con, w);                                  
        k3 = f(st + (par.s_step/st(1))/2*k2, con, w);         
        k4 = f(st + (par.s_step/st(1))*k3, con, w);           

        % x(k+1) = 
        % = f(x_opt(k),u_opt(k),w(k)) =  
        % = x_opt(k) + dt/6 * (k1(x_opt(k),u_opt(k)) + 2*k2(x_opt(k),u_opt(k)) + 2*k3(x_opt(k),u_opt(k)) + k4(x_opt(k),u_opt(k)))
        st_next_RK4= st + (par.s_step/st(1))/6*(k1 +2*k2 +2*k3 +k4);                  
        
        % variable to optimize x_opt(k+1), 
        st_next = X(:,k+1);                                                   
        
        % constraint x_opt(k+1) = x_opt(k) + dt/6 * (k1 +2*k2 +2*k3 +k4);                  
        % equal to x_opt(k+1) - x_opt(k) + dt/6 * (k1 +2*k2 +2*k3 +k4) == 0
        g.st_next = [g.st_next; st_next-st_next_RK4];                        
    end


    % compose final constraints vector
    g_nlp = [g_nlp; g.x0; g.st_next];
  
    %% Set boundary constraints to constraints vector
    % initialize "args" struct, to be used in the solver
    args = struct;
    
    % Initial condition constraint
    % x_opt(0) - x(0) == 0
    % Multiple shooting constraints
    % x_opt(k+1) - x_opt(k) + dt/6 * (k1 +2*k2 +2*k3 +k4) == 0
    args.lbg(1:par.n_states*(par.N+1)) = 0;                                                 
    args.ubg(1:par.n_states*(par.N+1)) = 0;  
    
    %% initialize boundary constraints to optimization variables vector

    % state velocity
    args.lbx(1:par.n_states:par.n_states*(par.N+1),1) = 0;                    
    args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = 40;                     
    
    % state battery energy
    args.lbx(2:par.n_states:par.n_states*(par.N+1),1) = 0.10 *par.E_bat_max;  
    args.ubx(2:par.n_states:par.n_states*(par.N+1),1) = 1.001 * par.E_bat_max;      
    
    % state time
    args.lbx(3:par.n_states:par.n_states*(par.N+1),1) = -Inf;                      
    args.ubx(3:par.n_states:par.n_states*(par.N+1),1) = Inf;                    
    
    % control input electric motor
    args.lbx(par.n_states*(par.N+1)+1:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = -5000; 
    args.ubx(par.n_states*(par.N+1)+1:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = 5000; 
    
    % control input brake (no regen)
    args.lbx(par.n_states*(par.N+1)+2:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = -5000;
    args.ubx(par.n_states*(par.N+1)+2:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = 0; 
    
    % slack variable S1, > 0, not used in the nlp
    args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = 0; 
    args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = 0; 
    
    % S2
    args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = 0; 
    args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = 0; 

end
