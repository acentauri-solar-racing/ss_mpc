syms P_t v_now v_past s m rho A_f c_d c_r g sin_alpha cos_alpha

F_r = c_r*m*g*cos_alpha;
F_g = m*g*sin_alpha;

v_mean = (v_now + v_past)/2;
t = s/v_mean;
a_mean = (v_now - v_past)/t;
eqn = P_t/v_mean - m*a_mean - 0.5*rho*A_f*c_d*v_mean^2 - F_r - F_g == 0;

v = solve(eqn, v_now, 'Real', true)
