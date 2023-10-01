%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 06.09.2023
%% Description 
% This function initializes all the variables needed for the mpc simulation
% loop, run the actual simulation and stores the results 

% INPUT: 
% "par": parameters struct,
% "args": argument struct needed for solver
% "f": symbolic function of states dynamics used in function "shift"
% "solver": solver to find optimal states and control inputs

% OUTPUT : 
% "par": parameters struct, added new parameters
% "OptResNLP": Optimization results
% "OptResNLP.xx": states at every step
% "OptResNLP.u_cl": control inputs at every step
% "OptResNLP.xx1": computed states trajectory

% "OptResNLP.xS": slack variable at every step 

% "OptResNLP.xG": irradiation at every step
% "OptResNLP.xwind_front": front wind at every step
% "OptResNLP.xwind_side": side wind at every step
% "OptResNLP.xtheta": ambient temperature at every step
% "OptResNLP.xSoC_N": Battery energy at horizon N, at every step
% "OptResNLP.xSoC_diff": Difference between Battery Energy(k+N) and Battery Energy target(k+N)
% "OptResNLP.xdist": distance vector

function [par, args, OptResNLP] = run_simulation_nlp(par, args, f, solver, s_0, t_0)
    % Start MPC
    %% NLP optimization over entire horizon
    par.N_NLP_t = par.N_NLP*par.s_step/15   % [s] Horizon length in seconds

    %% Initial conditions
    % overall iteration position, s = 0 => iter_initial = 0
    par.iter_initial = round(s_0/par.s_step);
    
    % initial velocity
    v_0 = par.Route.max_v(par.iter_initial+1)*0.95;
    
    % initial state of charge (value took from SoC target)
    SoC_0 = par.E_bat_target_DP(1+par.iter_initial);
    
    % initial conditions vector
    x0 = [v_0 ; SoC_0 ; t_0];  
    
    %initial slack value
    S0 =  0;                                
    
    % initialize optimal input vector
    OptResNLP.u_cl=[];     
    
    % initialize state save vector and slack save vector
    OptResNLP.xx(:,1) = x0;                    
    OptResNLP.xS(1) = S0;                      
    
    % initialize state trajectory vector for every MPC iteration
    OptResNLP.xx1 = zeros(par.N_NLP+1,par.n_states);                       
    OptResNLP.xdist(1) = s_0;
    
    % initialize trajectory control input vector (used as a "warm start" for
    % solver)
    u0 = zeros(par.N_NLP,par.n_controls);     
    
    % initialize trajectory states vector (used as a "warm start" for
    % solver)
    X0 = repmat(x0,1,par.N_NLP+1)';   

    %% Initialize Weather Data
    
    % initialize road inclination vector
    simvar.alpha = par.Route.incl';                                   
    
    % initialize battery energy target vector
    SoC_target = par.E_bat_target_DP;
    
    % initialize polynomial fit of G, fW
    [par.G_1, par.G_2, par.G_3] = get_poly(par.G_nlp, 0,  0+60*60*4);
    [par.fW_1, par.fW_2, par.fW_3] = get_poly(par.fW_nlp, 0,  0+60*60*4);
    [par.sW_1, par.sW_2, par.sW_3] = get_poly(par.sW_nlp, 0,  0+60*60*4);
    [par.temp_1, par.temp_2, par.temp_3] = get_poly(par.temp_nlp, 0,  0+60*60*4);

    % initialize parameters/prediction (warm start)
    vars_update_pred = [simvar.alpha(par.iter_initial+1:par.iter_initial+1+(par.N_NLP-1)); 
                           par.G_1;
                           par.G_2;
                           par.G_3;
                           par.fW_1;
                           par.fW_2;
                           par.fW_3;
                           par.sW_1;
                           par.sW_2;
                           par.sW_3;
                           par.temp_1;
                           par.temp_2;
                           par.temp_3
                           ];
    
    % initialize minimal/maximal velocity constraint
    args.lbx(1:par.n_states:par.n_states*(par.N_NLP+1),1) = 55/3.6; %par.Route.min_v(par.iter_initial+1:par.iter_initial+par.N_NLP+1);
    args.ubx(1:par.n_states:par.n_states*(par.N_NLP+1),1) = par.Route.max_v(par.iter_initial+1:par.iter_initial+par.N_NLP+1);
    
    %% Initialize Randbedingungen
    % state velocity
    v_idx = 1:par.n_states:par.n_states*(par.N_NLP+1);
    E_bat_idx = 2:par.n_states:par.n_states*(par.N_NLP+1);

    args.lbx(v_idx(end),1) = par.final_velocity_mpc-0.1;                  
    args.ubx(v_idx(end),1) = par.final_velocity_mpc+0.1;                     
    
    % state battery energy
    args.lbx(E_bat_idx(end),1) = par.final_E_bat_mpc-par.E_bat_max*0.00001;  
    args.ubx(E_bat_idx(end),1) = par.final_E_bat_mpc+par.E_bat_max*0.00001;     

    % slack variable, > 0
%     args.lbx(end,1) = 0; 
%     args.ubx(end,1) = 0; 
    %% Simulation Loop
    
    main_loop = tic;
    
    % update actual position 
    OptResNLP.xdist(1) = s_0;
    
    % update parameters "P" vector
    args.p = [x0; vars_update_pred; SoC_target(par.iter_initial+par.N_NLP)]; 
   
    % update optimization variables vector
    args.x0  = [reshape(X0',par.n_states*(par.N_NLP+1),1);reshape(u0',par.n_controls*par.N_NLP,1); S0]; 

    % run solver optimization
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    
    % get optimal input trajectory
    u = reshape(full(sol.x(par.n_states*(par.N_NLP+1)+1:end-1))',par.n_controls,par.N_NLP)';
    % store first optimal input
    OptResNLP.u_cl= [u];
    
    % store optimal states trajectory
    OptResNLP.xx1(:,1:3)= reshape(full(sol.x(1:par.n_states*(par.N_NLP+1)))',par.n_states,par.N_NLP+1)'; 
    
    % get slack variable value
    OptResNLP.S0 = reshape(full(sol.x(end))',1,1)';
    
    % store battery energy target at horizon N
    OptResNLP.xSoC_N(1) = SoC_target(par.iter_initial+par.N_NLP);
    
    main_loop_time = toc(main_loop);
    OptResNLP.average_mpc_time = main_loop_time
end

