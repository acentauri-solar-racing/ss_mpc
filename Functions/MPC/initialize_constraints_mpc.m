%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 06.09.2023
%% Description 
% This function initializes the constraints used in the MPC.
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
% "S": symbolic initialization for slack variable (used for loosening SoC
% target constraint),

% OUTPUT : 
% "par": new parameters added,
% "g_nlp": constraints vector,
%          "g,x0": constraint that enforces initial conditions,
%          "g.st_next": constraints of multiple shooting, where it is
%          enforced x(k+1) = f(x(k),u(k),w(k))
%          N.B. x(k), k in [1, N+1] are optimization variables (lifting)
%          "g.SoC_target": constraint that enforces the predicted SoC to be
%          greater or equal to the target SoC AT THE END of the prediction
%          horizon N; S is slack variable, S != 0 means that the constraint
%          is softened
% "args": struct used in optimizer, variable that can be updated during
%         simulation,
%         "args.lbg": contains lower bounds for constraints vector "g_nlp"
%         "args.ubg": contains upper bounds for constraint vector "g_nlp"
%         "args.lbx": contains lower bounds for optimization variables vector 
%         "args.ubx": contains upper bounds for optimization variables vector 
%          N.B. optimization variables vector contains states prediction,
%          control input prediction, slack variable

%%
function [par, g_nlp, args] = initialize_constraints_mpc(par, f, X, U, P, S1, S2)
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
    
    % constraint State of Charge target at end of Horizon N 
    % SoC_opt(k+N) - SoC_target(k+N) > 0, hard constraint
    % SoC_opt(k+N) - SoC_target(k+N) + S > 0, S >= 0, soft constraint
    g.SoC_target = [X(2,par.N+1) - P(par.n_states + par.N + 12 +1) + S1];     
    
    g.max_v = [];
    for k = 1:par.N+1
        g.max_v = [g.max_v; X(1,k) - S2];
    end

    % compose final constraints vector
    g_nlp = [g_nlp; g.x0; g.st_next; g.SoC_target; g.max_v];
  
    %% Set boundary constraints to constraints vector
    % initialize "args" struct, to be used in the solver
    args = struct;
    
    % Initial condition constraint
    % x_opt(0) - x(0) == 0
    % Multiple shooting constraints
    % x_opt(k+1) - x_opt(k) + dt/6 * (k1 +2*k2 +2*k3 +k4) == 0
    args.lbg(1:par.n_states*(par.N+1)) = 0;                                                 
    args.ubg(1:par.n_states*(par.N+1)) = 0;  

    %SoC has to be greater or equal than SoC target at the end of Horizon
    % SoC(k+N) > SoC_target(k+N), 
    args.lbg(par.n_states*(par.N+1)+1) =   0;      
    args.ubg(par.n_states*(par.N+1)+1) =   inf;

    args.lbg(par.n_states*(par.N+1)+2:par.n_states*(par.N+1)+1+par.N+1) = 0;                                                 
    args.ubg(par.n_states*(par.N+1)+2:par.n_states*(par.N+1)+1+par.N+1) = inf;  
    
    %% initialize boundary constraints to optimization variables vector
    % constraints values are initialized in the get_car_param file in
    % ss_offline_data

    % state velocity constraints (they will be updated in "run_simulation")
    args.lbx(1:par.n_states:par.n_states*(par.N+1),1) = 0;                    
    args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = 40;                     
    
    % state battery energy
    args.lbx(2:par.n_states:par.n_states*(par.N+1),1) = par.SoC_min *par.E_bat_max;  
    args.ubx(2:par.n_states:par.n_states*(par.N+1),1) = par.SoC_max * par.E_bat_max;      
    
    % state time
    args.lbx(3:par.n_states:par.n_states*(par.N+1),1) = -Inf;                      
    args.ubx(3:par.n_states:par.n_states*(par.N+1),1) = Inf;                    
    
    % control input electric motor
    args.lbx(par.n_states*(par.N+1)+1:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = par.P_el_min; 
    args.ubx(par.n_states*(par.N+1)+1:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = par.P_el_max; 
    
    % control input brake (no regen)
    args.lbx(par.n_states*(par.N+1)+2:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = par.P_brake_min;
    args.ubx(par.n_states*(par.N+1)+2:par.n_controls:par.n_states*(par.N+1)+par.n_controls*par.N,1) = par.P_brake_max;    
    
    % slack variable S1, > 0
    args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = 0; 
    args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = inf; 
    % S2
    args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = 0; 
    args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = inf; 

end
