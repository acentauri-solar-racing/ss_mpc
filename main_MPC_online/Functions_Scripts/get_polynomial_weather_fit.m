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
function [p_G, p_fW, p_rho, p_temp] = get_polynomial_weather_fit(weather, par, t_0, N_t)
    t_real = t_0;
    t_f_real = t_0 + N_t;
    t_initial_data = weather.timeline(1);
    t_0_data = t_real - t_initial_data;
    t_f_data = t_0_data + N_t;

    if t_0_data < 0
        t_0_data = 0;
    end

    G = weather.G_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),1);
    fW = weather.fW_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),1);
    rho = weather.rho_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),1);
    temp = weather.temp_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),1);
    t_x = linspace(floor(t_real/60/15),ceil(t_f_real/60/15),length(G));

    p_G = polyfit(t_x,G,2);
    p_fW = polyfit(t_x,fW,2);
    p_rho = polyfit(t_x,rho,2);
    p_temp = polyfit(t_x,temp,2);

end