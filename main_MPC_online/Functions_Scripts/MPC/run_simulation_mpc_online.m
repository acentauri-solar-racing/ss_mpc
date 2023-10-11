%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% This function initializes all the variables needed for the mpc simulation
% loop, run the actual simulation and stores the results 

% INPUT: 
% "par" (struct): parameters struct,
% "args" (struct): argument struct needed for solver
% "f" (function): symbolic function of states dynamics used in function "shift"
% "solver": solver to find optimal states and control inputs

% OUTPUT : 
% "par" (struct): parameters struct, added new parameters
% "OptRes" (struct): Optimization results

% "OptRes.xx": states at every step
% "OptRes.u_cl": control inputs at every step
% "OptRes.xx1": computed states trajectory

% "OptRes.xS1": slack variable at every step 

% "OptRes.xG(r)": irradiation at every step predicted by the controller (r instead is the "real" update) 
% "OptRes.xwind_front(r)": front wind at every step
% "OptRes.xwind_side(r)": side wind at every step
% "OptRes.xtheta(r)": ambient temperature at every step
% "OptRes.xSoC_N(r)": Battery energy at horizon N, at every step
% "OptRes.xSoC_diff(r)": Difference between Battery Energy(k+N) and Battery Energy target(k+N)
% "OptRes.xdist(r)": distance vector

function [par, OptRes] = run_simulation_mpc_online(par, weather, args, f, solver, s_0, DP_s_0, t_0)
    % Start MPC
    
    % overall iteration position, s = 0 => iter_initial = 0
    par.iter_initial = round(s_0/par.s_step);
    
    % mpc iteration
    par.iter_mpc = 0;             
    
    % maximum number of iterations
    par.iter_mpc_max = par.s_tot/par.s_step;
    
    % initial velocity
    v_0 = par.route.max_v(par.iter_initial+par.iter_mpc+1)*0.95;
    
    % initial state of charge (value took from SoC target)
    SoC_0 = par.E_bat_target_DP(1+par.iter_initial);
    
    % initial conditions vector
    x0 = [v_0 ; SoC_0 ; t_0];  
    
    %initial slack value
    S1_0 =  0;                                
    S2_0 = 0;
    % initialize optimal input vector
    OptRes.u_cl=[];     
    
    % initialize state save vector and slack save vector
    OptRes.xx(:,1) = x0;                    
    OptRes.xS1(1) = S1_0;   
    OptRes.xS2(1) = S2_0;                      

    
    % initialize state trajectory vector for every MPC iteration
    OptRes.xx1 = zeros(par.N+1,par.n_states, par.s_tot/par.s_step);                       
    OptRes.xdist(1) = s_0;
    
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
    [p_G, p_fW, p_rho, p_temp] = get_polynomial_weather_fit(weather, par, t_0, par.N_t);
    [p_G_r, p_fW_r, p_rho_r, p_temp_r] = get_polynomial_weather_fit(weather, par, t_0, 0);


    % initialize parameters/prediction (warm start)
    vars_update_pred = [par.route.incl(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.iter_mpc+1+(par.N-1)); 
                           p_G(1);
                           p_G(2);
                           p_G(3);
                           p_fW(1);
                           p_fW(2);
                           p_fW(3);
                           p_rho(1);
                           p_rho(2);
                           p_rho(3);
                           p_temp(1);
                           p_temp(2);
                           p_temp(3)
                           ];
    
    % initialize minimal/maximal velocity constraint
    args.lbx(1:par.n_states:par.n_states*(par.N+1),1) = 55/3.6; %par.route.min_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc);
    args.ubg(par.n_states*(par.N+1)+2:par.n_states*(par.N+1)+1+par.N+1) = par.route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc);
    args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)+10;                     

    %% Simulation Loop
    OptRes.iter_time = [];

    iter_time_1 = 0;
    iter_time_2 = 0;
    OptRes.skip = 0;

    main_loop = tic;
    change = 1;

    while par.iter_mpc < par.iter_mpc_max
        iter_time_1 = tic;
        iter_time_3 = tic;
        % update actual position 
        OptRes.xdist(par.iter_mpc+1) = s_0+par.s_step;
        
        % update parameters "P" vector
        args.p = [x0; vars_update_pred; par.SoC_target(par.N+par.iter_initial+par.iter_mpc)]; 
       
        % update optimization variables vector
        args.x0  = [reshape(X0',par.n_states*(par.N+1),1);reshape(u0',par.n_controls*par.N,1); S1_0; S2_0]; 
    
        % run solver optimization
        diary(['diagnostic/diary_' num2str(par.iter_mpc+1)])
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        diary off
        
        %%
        % get optimal input trajectory
        u = reshape(full(sol.x(par.n_states*(par.N+1)+1:end-2))',par.n_controls,par.N)';
        % store first optimal input
        %OptRes.ux1(:,1:2,par.iter_mpc+1) = reshape(full(sol.x(par.n_states*(par.N+1)+1:end-2))',par.n_controls,par.N)';
        OptRes.u_trajectory = u;
        OptRes.u_cl= [OptRes.u_cl ; u(1,:)];
        
        % store optimal states trajectory
        OptRes.xx1(:,1:3,par.iter_mpc+1)= reshape(full(sol.x(1:par.n_states*(par.N+1)))',par.n_states,par.N+1)'; 
        OptRes.trajectory = reshape(full(sol.x(1:par.n_states*(par.N+1)))',par.n_states,par.N+1)';

        % store slack variable value
        OptRes.xS1(par.iter_mpc+1) = reshape(full(sol.x(end-1))',1,1)';
        OptRes.xS2(par.iter_mpc+1) = reshape(full(sol.x(end))',1,1)';
        
        %%
        % store weather condition Data (perfect fit)
        OptRes.xGr(par.iter_mpc+1) = p_G_r(1)*(x0(3)/60/15)^2 + p_G_r(2)*(x0(3)/60/15) + p_G_r(3);
        OptRes.xfWr(par.iter_mpc+1) = p_fW_r(1)*(x0(3)/60/15)^2 + p_fW_r(2)*(x0(3)/60/15) + p_fW_r(3);
        OptRes.xrhor(par.iter_mpc+1) = p_rho_r(1)*(x0(3)/60/15)^2 + p_rho_r(2)*(x0(3)/60/15) + p_rho_r(3);
        OptRes.xtempr(par.iter_mpc+1) = p_temp_r(1)*(x0(3)/60/15)^2 + p_temp_r(2)*(x0(3)/60/15) + p_temp_r(3);

        % store weather condition Data (fit for mpc)
        OptRes.xG(par.iter_mpc+1) = p_G(1)*(x0(3)/60/15)^2 + p_G(2)*(x0(3)/60/15) + p_G(3);
        OptRes.xfW(par.iter_mpc+1) = p_fW(1)*(x0(3)/60/15)^2 + p_fW(2)*(x0(3)/60/15) + p_fW(3);
        OptRes.xrho(par.iter_mpc+1) = p_rho(1)*(x0(3)/60/15)^2 + p_rho(2)*(x0(3)/60/15) + p_rho(3);
        OptRes.xtemp(par.iter_mpc+1) = p_temp(1)*(x0(3)/60/15)^2 + p_temp(2)*(x0(3)/60/15) + p_temp(3);
        
        %%
        % store battery energy target at horizon N
        %OptRes.xSoC_N(par.iter_mpc+1) = par.SoC_target(par.iter_initial+par.N+par.iter_mpc);
    
        % store battery energy difference between prediction and target
        OptRes.xSoC_diff(par.iter_mpc+1) = OptRes.xx1(end,2,par.iter_mpc+1) - par.SoC_target(par.iter_initial+par.N+par.iter_mpc);
        
        iter_time_1 = toc(iter_time_1);
        %%

        
        if iter_time_1+iter_time_2-0.1 < par.s_step/x0(1)
            pause(par.s_step/x0(1)-iter_time_1-iter_time_2)
        else
            disp('mpc slower than car')
            OptRes.skip = OptRes.skip+1;
        end
        
        iter_time_2 = tic;
        
        % advance simulation, update real plant
        [s_0, x0, u0] = shift(par.s_step, s_0, x0, u, f, [par.route.incl(par.iter_initial+par.iter_mpc+1); ...
                                                           p_G_r(1);
                                                           p_G_r(2);
                                                           p_G_r(3);
                                                           p_fW_r(1);
                                                           p_fW_r(2);
                                                           p_fW_r(3);
                                                           p_rho_r(1);
                                                           p_rho_r(2);
                                                           p_rho_r(3);
                                                           p_temp_r(1);
                                                           p_temp_r(2);
                                                           p_temp_r(3)]);
        %%
    
        % give real online update

        %%
        [p_G, p_fW, p_rho, p_temp] = get_polynomial_weather_fit(weather, par, x0(3), par.N_t);
        [p_G_r, p_fW_r, p_rho_r, p_temp_r] = get_polynomial_weather_fit(weather, par, x0(3), 0);

        % update the weather/road variables
        vars_update_pred = [par.route.incl(par.iter_initial+par.iter_mpc+1 +1:par.iter_initial+par.iter_mpc+1+(par.N-1) +1); 
                            p_G(1);
                            p_G(2);
                            p_G(3);
                            p_fW(1);
                            p_fW(2);
                            p_fW(3);
                            p_rho(1);
                            p_rho(2);
                            p_rho(3);
                            p_temp(1);
                            p_temp(2);
                            p_temp(3)];     
        
        %%
        % update the maximal velocity constraint
        args.ubg(par.n_states*(par.N+1)+2:par.n_states*(par.N+1)+1+par.N+1) = par.route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc);
        args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)+10;                     
        
        %% ADD disturbances
        % x0(1) = x0(1) + rand()*1/3.6 - rand()*1/3.6;
        % x0(2) = x0(2) + rand()*par.E_bat_max*0.000001 - rand()*par.E_bat_max*0.000001;

%         if par.iter_mpc == 30;
%             x0(1) = 60/16.6;
%             x0(2) = par.E_bat_max*0.99;
%         end


        %% store real state values
        OptRes.xx(:,par.iter_mpc+2) = x0;                     
        
        %% Online Plot
        plotOptTrajectory(OptRes, par, s_0, x0)
        
        %%
        % get states trajectory to use as "args.x0"
        X0 = reshape(full(sol.x(1:par.n_states*(par.N+1)))',3,par.N+1)'; 
       
        % Shift trajectory to initialize the next step
        X0 = [X0(2:end,:);X0(end,:)];

        % update iter_mpc
        par.iter_mpc
        par.iter_mpc = par.iter_mpc + 1;
        iter_time_2 = toc(iter_time_2);
        iter_time_3 = toc(iter_time_3);
        OptRes.iter_time = [OptRes.iter_time; iter_time_3];

        if kbhit()
            break;
        end
    end
    
    main_loop_time = toc(main_loop);
    OptRes.online_time = main_loop_time
    OptRes.average_mpc_time = main_loop_time/(par.iter_mpc+1)

end