%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 08.09.2023
%% Description 
% This function returns the polinomial fit parameters for the solar
% irradiation G, assuming it's only depending on the time

% INPUT: 
% "par": parameters struct
% "t_0": initial time of polinomial fit
% "t_f": final time of polinomial fit


% OUTPUT : 
% G_1, G_2, G_3: second order polinomial fit parameters, G_1*t^2 + G_2*t + G_3

%%
function [p_G, p_fW, p_rho, p_temp] = get_polynomial_weather_fit(weather, par, N, t_0, N_t)
    
    %real time, machine time
    t_real = t_0;
    % final time prediction
    t_f_real = t_0 + N_t;
    % initial time of data
    t_initial_data = weather.timeline(1);
    % difference between real time and initial time of data, to understand
    % which initial point of data to use
    t_0_data = t_real - t_initial_data;
    % final time of the data
    t_f_data = t_0_data + N_t;

    if t_0_data < 0
        t_0_data = 0;
    end
    
    % Find the index of the nearest distance between cumulative distance,
    % s_0 and final position s_0+s_step*N
    [~, idx_0] = min(abs(weather.cumdist - par.s_0));
    [~, idx_end] = min(abs(weather.cumdist - par.s_0+N*par.s_step));
    
    % Cut data from t_0 to t_f, take spatial row mean from s_0 to
    % s_0+s_step*N
    G = weather.G_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_0:idx_end);
    G_mean = mean(G,2);
    fW = weather.fW_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_0:idx_end);
    fW_mean = mean(fW,2);
    rho = weather.rho_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_0:idx_end);
    rho_mean = mean(rho,2);
    temp = weather.temp_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_0:idx_end);
    temp_mean = mean(temp,2);
    t_x = linspace(floor(t_real/60/15),ceil(t_f_real/60/15),length(G));

    p_G = polyfit(t_x,G_mean,2);
    p_fW = polyfit(t_x,fW_mean,2);
    p_rho = polyfit(t_x,rho_mean,2);
    p_temp = polyfit(t_x,temp_mean,2);

end