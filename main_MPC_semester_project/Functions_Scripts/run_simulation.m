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
% "OptRes": Optimization results
% "OptRes.xx": states at every step
% "OptRes.u_cl": control inputs at every step
% "OptRes.xx1": computed states trajectory

% "OptRes.xS": slack variable at every step 

% "OptRes.xG": irradiation at every step
% "OptRes.xwind_front": front wind at every step
% "OptRes.xwind_side": side wind at every step
% "OptRes.xtheta": ambient temperature at every step
% "OptRes.xSoC_N": Battery energy at horizon N, at every step
% "OptRes.xSoC_diff": Difference between Battery Energy(k+N) and Battery Energy target(k+N)
% "OptRes.xdist": distance vector

function [par, OptRes] = run_simulation(par, args, f, solver, s_0, t_0)
    % Start MPC
    
    % overall iteration position, s = 0 => iter_initial = 0
    par.iter_initial = round(s_0/par.s_step);
    
    % mpc iteration
    par.iter_mpc = 0;             
    
    % maximum number of iterations
    par.iter_mpc_max = par.s_tot/par.s_step;
    
    % initial velocity
    v_0 = par.Route.max_v(par.iter_initial+par.iter_mpc+1)*0.95;
    
    % initial state of charge (value took from SoC target)
    SoC_0 = par.E_bat_target_DP(1+par.iter_initial);
    
    % initial conditions vector
    x0 = [v_0 ; SoC_0 ; t_0];  
    
    %initial slack value
    S0 =  0;                                
    
    % initialize optimal input vector
    OptRes.u_cl=[];     
    
    % initialize state save vector and slack save vector
    OptRes.xx(:,1) = x0;                    
    OptRes.xS(1) = S0;                      
    
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
    simvar.alpha = par.Route.incl';                                   
    
    % initialize battery energy target vector
    SoC_target = par.E_bat_target_DP;
    
    % initialize polynomial fit of G used for the MPC prediction
    [par.G_1, par.G_2, par.G_3] = get_G(par.G_MPC, 0, 0+ par.N*par.s_step/33.3);

    % initialize polynomial fit of G used for real conditions shift
    [par.G_1_shift, par.G_2_shift, par.G_3_shift] = get_G(par.G_real, 0, 0+ 180);

    % initialize parameters/prediction (warm start)
    vars_update_pred = [simvar.alpha(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.iter_mpc+1+(par.N-1)); 
                        linspace(0,0, par.N)';                    
                        linspace(0,0, par.N)'; 
                        linspace(25, 25, par.N)';
                        par.G_1;
                        par.G_2;
                        par.G_3];
    
    % initialize minimal/maximal velocity constraint
    args.lbx(1:par.n_states:par.n_states*(par.N+1),1) = 55/3.6; %par.Route.min_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc);
    args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc);
    
    %% Simulation Loop
    
    main_loop = tic;
    while par.iter_mpc < par.iter_mpc_max
    
        % update actual position 
        OptRes.xdist(par.iter_mpc+1) = s_0;
        
        % update parameters "P" vector
        args.p = [x0; vars_update_pred; SoC_target(par.iter_initial+par.N+par.iter_mpc)]; 
       
        % update optimization variables vector
        args.x0  = [reshape(X0',par.n_states*(par.N+1),1);reshape(u0',par.n_controls*par.N,1); S0]; 
    
        % run solver optimization
        diary(['diagnostic/diary_' num2str(par.iter_mpc+1)])
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        diary off
        
        % get optimal input trajectory
        u = reshape(full(sol.x(par.n_states*(par.N+1)+1:end-1))',par.n_controls,par.N)';
        % store first optimal input
        OptRes.u_cl= [OptRes.u_cl ; u(1,:)];
        
        % store optimal states trajectory
        OptRes.xx1(:,1:3,par.iter_mpc+1)= reshape(full(sol.x(1:par.n_states*(par.N+1)))',par.n_states,par.N+1)'; 
        
        % get slack variable value
        S0 = reshape(full(sol.x(end))',1,1)';
        % store slack variable value
        OptRes.xS(par.iter_mpc+1) = S0;
        
        % store weather condition Data
        OptRes.xG_MPC(par.iter_mpc+1) = par.G_1*(x0(3)/60)^2 + par.G_2*(x0(3)/60) + par.G_3;
        OptRes.xG(par.iter_mpc+1) = par.G_1_shift*(x0(3)/60)^2 + par.G_2_shift*(x0(3)/60) + par.G_3_shift;
        OptRes.xwind_front(par.iter_mpc+1) = par.wind_front_int({s_0,x0(3)/60/60});
        OptRes.xwind_side(par.iter_mpc+1) = par.wind_side_int({s_0,x0(3)/60/60});
        OptRes.xtheta(par.iter_mpc+1) = par.theta_int({s_0,x0(3)/60/60});
        
        % store battery energy target at horizon N
        OptRes.xSoC_N(par.iter_mpc+1) = SoC_target(par.iter_initial+par.N+par.iter_mpc);
    
        % store battery energy difference between prediction and target
        OptRes.xSoC_diff(par.iter_mpc+1) = OptRes.xx1(end,2,par.iter_mpc+1) - SoC_target(par.iter_initial+par.N+par.iter_mpc);
        

        % advance simulation, update real plant
        [s_0, x0, u0] = shift(par.s_step, s_0, x0, u, f, [simvar.alpha(par.iter_initial+par.iter_mpc+1); ...          
                                           0; ...
                                           0; ...
                                           25; ...
                                           par.G_1_shift; ...
                                           par.G_2_shift; ...
                                           par.G_3_shift]);
    
        % update G polynomial fit
        [par.G_1, par.G_2, par.G_3] = get_G(par.G_MPC, x0(3), x0(3) + 600);
        % update real G polynomial
        [par.G_1_shift, par.G_2_shift, par.G_3_shift] = get_G(par.G_real, x0(3), x0(3)+ 180);

        % update the weather/road variables for the MPC
        vars_update_pred = [simvar.alpha(par.iter_initial+par.iter_mpc+1 +1:par.iter_initial+par.iter_mpc+1+(par.N-1) +1); 
                           linspace(0,0, par.N)';                    
                           linspace(0,0, par.N)'; 
                           linspace(25, 25, par.N)';
                           par.G_1;
                           par.G_2;
                           par.G_3];     
    
        % update the maximal velocity constraint
        %args.lbx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.min_v(par.iter_initial+par.iter_mpc+1 +1:par.iter_initial+par.N+1+par.iter_mpc +1); %state velocity lower bound
        args.ubx(1:par.n_states:par.n_states*(par.N+1),1) = par.Route.max_v(par.iter_initial+par.iter_mpc+1 +1:par.iter_initial+par.N+1+par.iter_mpc +1); %state velocity upper bound
    
    
        % store real state values
        OptRes.xx(:,par.iter_mpc+2) = x0;                     
        % get states trajectory to use as "args.x0"
        X0 = reshape(full(sol.x(1:par.n_states*(par.N+1)))',3,par.N+1)'; 
       
        % Shift trajectory to initialize the next step
        X0 = [X0(2:end,:);X0(end,:)];

        % update iter_mpc
        par.iter_mpc
        par.iter_mpc = par.iter_mpc + 1;
    end
    
    main_loop_time = toc(main_loop);
    OptRes.average_mpc_time = main_loop_time/(par.iter_mpc+1)
end