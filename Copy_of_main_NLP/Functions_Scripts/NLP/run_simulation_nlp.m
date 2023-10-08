%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% This function initializes all the variables needed for the nlp simulation
% run the actual simulation and stores the results 

% INPUT: 
% "par": parameters struct,
% "args": argument struct needed for solver
% "solver": solver to find optimal states and control inputs
% "s_0": in [m]; initial position 
% "t_0": in [s]; initial time (t = 0 should mean 08.00/start of the day
% "v_0": in [m/s]; initial velocity (t = 0 should mean 08.00/start of the day)
% "SoC_0": in [-]; initial SoC

% "final_velocity_MPC":in [m/s]; boundary condition for final velocity, uses MPC
% final value (to compare benchmark); use -1 if you wanna set an arbitrary
% value ([0.75-1]*route_max_velocity)
% "final_E_bat_MPC": in [J]; boundary condition for final battery energy, uses MPC
% final value (to compare benchmark); use -1 if you wanna set an arbitrary
% value (uses DP solution)

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

function [par, OptResNLP] = run_simulation_nlp(par, weather, args, solver, s_0, t_0, v_0, SoC_0, final_velocity_MPC, final_E_bat_MPC)
%% NLP optimization over entire horizon
    %% Initial conditions
    % overall iteration position, s = 0 => iter_initial = 0
    par.iter_initial = round(s_0/par.s_step);
    par.N_t = par.N*par.s_step/15; % prediction in time domain, how far the weather is interpolated, assuming minimum velocity is 30 km/h (8.33 m/s)

    
    if v_0 == -1
        % initial velocity
        v_0 = par.route.max_v(par.iter_initial+1)*0.95;
    end

    if SoC_0 == -1
        % initial state of charge (value took from SoC target)
        SoC_0 = par.E_bat_target_DP(1+par.iter_initial);
    end

    % initial conditions vector
    x0 = [v_0 ; SoC_0; t_0];  
    
    %initial slack value
    S1_0 =  0;                                
    S2_0 =  0;                                

    % initialize optimal input vector
    OptResNLP.u_cl=[];     
    
    % initialize state save vector and slack save vector
    OptResNLP.xx(:,1) = x0;                    
    OptResNLP.xS1(1) = S1_0;    
    OptResNLP.xS2(1) = S2_0;                      

    
    % initialize state trajectory vector for every MPC iteration
    OptResNLP.xx1 = zeros(par.N+1,par.n_states);                       
    OptResNLP.xdist(1) = s_0;
    
    % initialize trajectory control input vector (used as a "warm start" for
    % solver)
    u0 = zeros(par.N,par.n_controls);     
    
    % initialize trajectory states vector (used as a "warm start" for
    % solver)
    X0 = repmat(x0,1,par.N+1)';   

    %% Initialize Weather Data
    
    % initialize road inclination vector
    par.route.incl = par.route.incl';                                   
    
    % initialize battery energy target vector
    SoC_target = par.E_bat_target_DP;
    
    % initialize polynomial fit of G, fW
    [par.G_1, par.G_2, par.G_3, par.G_4] = get_poly(weather.G_data, 0, par.N_t);
    [par.fW_1, par.fW_2, par.fW_3, par.fW_4] = get_poly(weather.fW_data, 0, par.N_t);
    [par.sW_1, par.sW_2, par.sW_3, par.sW_4] = get_poly(weather.sW_data, 0, par.N_t);
    [par.temp_1, par.temp_2, par.temp_3, par.temp_4] = get_poly(weather.temp_data, 0, par.N_t);

    % initialize parameters/prediction (warm start)
    vars_update_pred = [par.route.incl(par.iter_initial+1:par.iter_initial+1+(par.N-1)); 
                           par.G_1;
                           par.G_2;
                           par.G_3;
                           par.G_4;
                           par.fW_1;
                           par.fW_2;
                           par.fW_3;
                           par.fW_4;
                           par.sW_1;
                           par.sW_2;
                           par.sW_3;
                           par.sW_4;
                           par.temp_1;
                           par.temp_2;
                           par.temp_3;
                           par.temp_4
                           ];
    
    % initialize minimal/maximal velocity constraint
    args.lbx(1:par.n_states:par.n_states*(par.N+1),1) = 50/3.6; %par.route.min_v(par.iter_initial+1:par.iter_initial+par.N+1);
    args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.route.max_v(par.iter_initial+1:par.iter_initial+par.N+1);
    
    %% Initialize Randbedingungen
    % state velocity
    v_idx = 1:par.n_states:par.n_states*(par.N+1);
    args.lbx(v_idx(end),1) = final_velocity_MPC -0.001;                    
    args.ubx(v_idx(end),1) = final_velocity_MPC +0.001;                     
    
    if final_velocity_MPC == -1
        args.lbx(v_idx(end),1) = par.route.max_v(par.iter_initial+par.N+1)*0.75;                    
        args.ubx(v_idx(end),1) = par.route.max_v(par.iter_initial+par.N+1); 
    end

    % state battery energy
    E_bat_idx = 2:par.n_states:par.n_states*(par.N+1);

    args.lbx(E_bat_idx(end),1) = final_E_bat_MPC -par.E_bat_max*0.00001;  
    % args.ubx(E_bat_idx(end),1) = final_E_bat_MPC +par.E_bat_max*0.00001;  

    if final_E_bat_MPC == -1
        args.lbx(E_bat_idx(end),1) = par.E_bat_target_DP(1+par.iter_initial+par.N) -par.E_bat_max*0.00001;  
        %args.ubx(E_bat_idx(end),1) = par.E_bat_target_DP(1+par.iter_initial+par.N) +par.E_bat_max*0.00001;  
        args.ubx(E_bat_idx(end),1) = par.E_bat_max;  
    end

    % slack variable, they are 0 in the NLP
    args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = 0; 
    args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = 0; 
    args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = 0; 
    args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = 0; 
    %% Simulation Loop
    
    % measure computational time
    main_loop = tic;
    
    % update actual position 
    OptResNLP.xdist(1) = s_0;
    
    % update parameters "P" vector
    args.p = [x0; vars_update_pred; SoC_target(par.iter_initial+par.N)]; 
   
    % update optimization variables vector
    args.x0  = [reshape(X0',par.n_states*(par.N+1),1);reshape(u0',par.n_controls*par.N,1); S1_0; S2_0]; 

    % run solver optimization
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    
    % get optimal input trajectory
    u = reshape(full(sol.x(par.n_states*(par.N+1)+1:end-2))',par.n_controls,par.N)';
    % store first optimal input
    OptResNLP.u_cl= [OptResNLP.u_cl ; u(:,:)];
    
    % store optimal states trajectory
    OptResNLP.xx1(:,1:3)= reshape(full(sol.x(1:par.n_states*(par.N+1)))',par.n_states,par.N+1)'; 
    
    % store battery energy target at horizon N
    OptResNLP.xSoC_N(1) = SoC_target(par.iter_initial+par.N);
    
    main_loop_time = toc(main_loop);
    OptResNLP.nlp_time = main_loop_time
end