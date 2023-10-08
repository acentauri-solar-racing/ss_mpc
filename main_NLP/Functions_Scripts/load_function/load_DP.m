function par = load_DP(par, filenamepath)
        par.E_bat_target_DP_raw = load(filenamepath).OptRes.states.E_bat*3600; % [Wh = W * 3600s = J]
        par.E_bat_target_DP = interp1(linspace(0,300,300+1), par.E_bat_target_DP_raw, linspace(0,300,(300)*round(par.s_step_DP/par.s_step)+1));
        
        par.v_DP_raw = load(filenamepath).OptRes.states.V; % [m/s]
        par.v_DP = interp1(linspace(0,300,300+1), par.v_DP_raw, linspace(0,300,(300)*round(par.s_step_DP/par.s_step)+1));
        
        par.P_mot_el_DP_raw = load(filenamepath).OptRes.inputs.P_mot_el; % [m/s]
        par.P_mot_el_DP = interp1(linspace(1,300,300), par.P_mot_el_DP_raw , linspace(1,300,300*round(par.s_step_DP/par.s_step)));
end