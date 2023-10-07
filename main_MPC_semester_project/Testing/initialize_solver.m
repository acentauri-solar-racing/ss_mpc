%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 06.09.2023
%% Description 
% This function creates the optimizer solver used in the simulation loop;
% the solver is IPOPT

% INPUT: 
% "par": parameters struct,
% "obj": cost function,
% *g_nlp": constraints vector,
% "X": symbolic initialization of states (velocity, battery energy ,time),
% "U": symbolic initialization of control input (Electric Motor, Brake),
% "P": symbolic initialization of parameters (v0, E_bat0, t0, road inclination,
% irradiation, front wind velocity, side wind velocity, ambient
% temperature, E_bat_target),
% "S": symbolic initialization for slack variable (used for loosening SoC
% target constraint),

% OUTPUT : 
% "solver": optimization solver to be used to find optimal states and
% controls

%
%%
function solver = initialize_solver(par, obj, g_nlp, X, U, P, S1, S2)
    %import casadi.* package
    import casadi.*
    
    % make the decision/optimization variable one column vector, to be
    % given at the solver
    OPT_variables = [reshape(X,par.n_states*(par.N+1),1);reshape(U,par.n_controls*par.N,1); S1; S2];           % OPT_variables = [x1, x2, x3, ... xN, u1, u2, u3, ..., uN-1]
    
    % create symbolic non linear problem
    nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g_nlp, 'p', P); 
    
    % solver options
    opts = struct;
    opts.ipopt.max_iter = 2000;
    opts.ipopt.print_level = 3;     %0,3,5
    opts.print_time = 0;
    opts.ipopt.acceptable_tol =1e-8; %1e-8
    opts.ipopt.acceptable_obj_change_tol = 1e-8;
    
    % create solver
    solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

end