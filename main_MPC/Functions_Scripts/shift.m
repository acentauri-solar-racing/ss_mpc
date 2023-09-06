%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 06.09.2023
%% Description 
% INPUT: 
% "s_step": discretization step
% "s_0": position at the moment before the shift
% "x0": initial conditions before the shift
% "u": control input
% "f": function of state dynamics
% "vars": varying condition (road inclination, irradiation, front wind,
% side wind, amb. temperature)
% OUTPUT : 
% "s_0": new position after shift
% "x0": new initial conditions after shift
% "u0": input prediction after shift

%%
function [s_0, x0, u0] = shift(s_step, s_0, x0, u, f, vars)
    st = x0;
    con = u(1,:)';
    
    k1 = f(st, con, vars); 
    k2 = f(st + (s_step/st(1))/2*k1, con, vars);     % reminder s_step = ds, st(1) = v, dt = ds/v
    k3 = f(st + (s_step/st(1))/2*k2, con, vars); 
    k4 = f(st + (s_step/st(1))*k3, con, vars); 
    st = st + (s_step/st(1))/6*(k1 +2*k2 +2*k3 +k4); %RK 4th order method    
    
    % f_value = f(st,con,vars);
    % st = st+ ((s_step/st(1))*f_value);               % dt = ds/dv;
    
    
    x0 = full(st);                              % update x0 with new state
    s_0 = s_0 + s_step;                                % update distance s_0 with new state
    u0 = [u(2:size(u,1),:);u(size(u,1),:)];         
end
