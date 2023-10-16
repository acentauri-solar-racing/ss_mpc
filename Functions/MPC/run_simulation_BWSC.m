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

% "velocity_target":in [m/s]; boundary condition for final velocity, uses MPC
% final value (to compare benchmark); use -1 if you wanna set an arbitrary
% value ([0.75-1]*route_max_velocity)
% "SoC_target": in [J]; boundary condition for final battery energy, uses MPC
% final value (to compare benchmark); use -1 if you wanna set an arbitrary
% value (uses DP solution)

% OUTPUT : 
% "par": parameters struct, added new parameters
% "OptRes": Optimization results
% "OptRes.u_cl": control inputs at every step
% "OptRes.xx": computed states trajectory

% "OptRes.xS": slack variable at every step 

% "OptRes.xG": irradiation at every step
% "OptRes.xwind_front": front wind at every step
% "OptRes.xwind_side": side wind at every step
% "OptRes.xtheta": ambient temperature at every step
% "OptRes.xSoC_N": Battery energy at horizon N, at every step
% "OptRes.xSoC_diff": Difference between Battery Energy(k+N) and Battery Energy target(k+N)
% "OptRes.xdist": distance vector

function [par, weather, OptRes] = run_simulation_BWSC(par, weather, args, solver, velocity_target, SoC_target)
%% NLP optimization over entire horizon
    %% Initial conditions
%% TO UPDATE NEW POSITION, VELOCITY AND SOC     
%     s_0 = get_initial_position();
%     v_0 = get_initial_velocity();
%     SoC_0 = get_initial_SoC();

    s_0 = par.s_0;
    t_0 = par.t_0;
    par.iter_initial = round(s_0/par.s_step);
    v_0 = par.route.max_v(par.iter_initial+1)*0.95;
    % initial state of charge (value took from SoC target)
    SoC_0 = par.E_bat_target_DP(1+par.iter_initial)/par.E_bat_max;

    % initial conditions vector
    x0 = [v_0 ; SoC_0*par.E_bat_max; t_0];  
%%

    %initial slack value
    S1_0 =  0;                                
    S2_0 =  0;                                

    % initialize optimal input vector
    OptRes.u_cl=[];     
    
    % initialize state save vector and slack save vector
    OptRes.xx(:,1) = x0;                    
%   OptRes.xS1(1) = S1_0;    
%   OptRes.xS2(1) = S2_0;                      

    
    % initialize state trajectory vector for every MPC iteration
    OptRes.xx = zeros(par.N+1,par.n_states);                       
    OptRes.xdist = linspace(s_0, s_0+(par.N)*par.s_step,par.N+1)';
    
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
    par.SoC_target = par.E_bat_target_DP;
    
    % initialize polynomial fit of G, fW
    [weather.p_G, weather.p_fW, weather.p_rho, weather.p_temp] = get_polynomial_weather_fit(weather, par, t_0, par.N_t);

    % initialize parameters/prediction (warm start)
    vars_update_pred = [par.route.incl(par.iter_initial+1:par.iter_initial+1+(par.N-1)); 
                           weather.p_G(1);
                           weather.p_G(2);
                           weather.p_G(3);
                           weather.p_fW(1);
                           weather.p_fW(2);
                           weather.p_fW(3);
                           weather.p_rho(1);
                           weather.p_rho(2);
                           weather.p_rho(3);
                           weather.p_temp(1);
                           weather.p_temp(2);
                           weather.p_temp(3)
                           ];
    
    % initialize minimal/maximal velocity constraint over N
    args.lbx(1:par.n_states:par.n_states*(par.N+1),1) = par.v_min; %par.route.min_v(par.iter_initial+1:par.iter_initial+par.N+1);
    args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.route.max_v(par.iter_initial+1:par.iter_initial+par.N+1);
    
    %% Initialize Randbedingungen
    % state velocity target at the end distance
    v_idx = 1:par.n_states:par.n_states*(par.N+1);
    args.lbx(v_idx(end),1) = velocity_target -0.001;                    
    args.ubx(v_idx(end),1) = velocity_target +0.001;                     
    
    if velocity_target == -1
        args.lbx(v_idx(end),1) = par.route.max_v(par.iter_initial+par.N+1)*par.v_end_lb_percentage;                    
        args.ubx(v_idx(end),1) = par.route.max_v(par.iter_initial+par.N+1)*par.v_end_ub_percentage; 
    end

    % state battery energy target at the end distance
    E_bat_idx = 2:par.n_states:par.n_states*(par.N+1);
    args.lbx(E_bat_idx(end),1) = SoC_target*par.E_bat_max -par.E_bat_max*0.00001;  

%     args.ubx(E_bat_idx(end),1) = par.E_bat_max;  
%     if SoC_target == -1
%         args.lbx(E_bat_idx(end),1) = par.E_bat_target_DP(1+par.iter_initial+par.N) -par.E_bat_max*0.00001;  
%         %args.ubx(E_bat_idx(end),1) = par.E_bat_target_DP(1+par.iter_initial+par.N) +par.E_bat_max*0.00001;  
%     end


    % slack variables
%     args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = 0; 
%     args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +1,1) = inf; 
%     args.lbx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = 0; 
%     args.ubx(par.n_states*(par.N+1)+par.n_controls*par.N +2,1) = 0; 
    %% Simulation Loop
    
    % measure computational time
    main_loop = tic;
    
    % update actual position 
    
    % update parameters "P" vector
    args.p = [x0; vars_update_pred; par.E_bat_target_DP(par.iter_initial+par.N)]; 
   
    % update optimization variables vector
    args.x0  = [reshape(X0',par.n_states*(par.N+1),1);reshape(u0',par.n_controls*par.N,1); S1_0; S2_0]; 

    % run solver optimization
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    
    % get optimal input trajectory
    u = reshape(full(sol.x(par.n_states*(par.N+1)+1:end-2))',par.n_controls,par.N)';
    slack = reshape(full(sol.x(end-1))',1,1)';
    % store first optimal input
    OptRes.u_cl= [OptRes.u_cl ; u(:,:)];
    
    % store optimal states trajectory
    OptRes.xx(:,1:3)= reshape(full(sol.x(1:par.n_states*(par.N+1)))',par.n_states,par.N+1)'; 
    
    % store battery energy target at horizon N
    OptRes.xSoC_N(1) = par.E_bat_target_DP(par.iter_initial+par.N);
    
    main_loop_time = toc(main_loop);
    OptRes.nlp_time = main_loop_time;
    
%     %% Print info on time elapsed in real life and simulation, computational time
      printresults.computational_time = toc(main_loop);
%     OptRes.online_time = main_loop_time;
% 
%     printresult.time_to_run_distance = (OptRes.xx(3,end)-par.t_0);              % [min]      args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)*1.1;                     
%     printresult.average_nlp_time_computation_iter = main_loop_time;
%     printresult.v_end = OptRes.xx(1,end);
%     printresult.SoC_end = OptRes.xx(2,end)/par.E_bat_max
    fprintf('slack variable = %d \n', slack);
    fprintf('computational time = %d s\n', main_loop_time);



end