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
    
    % Find the index of the nearest distance between cumulative distance
    % and middle point between s_0 final prediction position
    [~, idx_dist] = min(abs(weather.cumdist - par.s_0+par.N*par.s_step/2));
    
    G = weather.G_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_dist);
    fW = weather.fW_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_dist);
    rho = weather.rho_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_dist);
    temp = weather.temp_data(1+floor(t_0_data/60/15):1+ceil(t_f_data/60/15),idx_dist);
    t_x = linspace(floor(t_real/60/15),ceil(t_f_real/60/15),length(G));

    p_G = polyfit(t_x,G,2);
    p_fW = polyfit(t_x,fW,2);
    p_rho = polyfit(t_x,rho,2);
    p_temp = polyfit(t_x,temp,2);

    G_fit = polyval(p_G,t_x);
    fW_fit = polyval(p_fW,t_x);
    rho_fit = polyval(p_rho,t_x);
    temp_fit = polyval(p_temp,t_x);

    
    %%
    % Convert t_x to seconds
    t_x_seconds = t_x * 60 * 15;
    
    % Create a cell array of "HH:MM" formatted labels
    time_labels = cell(length(t_x), 1);
    for i = 1:length(t_x)
        hours = floor(t_x_seconds(i) / 3600);
        minutes = mod(floor(t_x_seconds(i) / 60), 60);
        time_labels{i} = sprintf('%02d:%02d', hours, minutes);
    end

    figure
    subplot(4,1,1)
    plot(t_x,G), hold on;
    plot(t_x,G_fit)
    title('Solar Irradiance G')
    legend('mean', 'mean fit', 'Location','northeast')
    % Set custom x-axis tick values and labels
    xticks(t_x);
    xticklabels(time_labels);

    subplot(4,1,2)
    plot(t_x, fW), hold on;
    plot(t_x, fW_fit)    
    title('Front Wind Velocity fW')
    xticks(t_x);
    xticklabels(time_labels);

    subplot(4,1,3)
    plot(t_x, rho), hold on;
    plot(t_x, rho_fit)    
    title('Air Density \rho')
    xticks(t_x);
    xticklabels(time_labels);

    subplot(4,1,4)
    plot(t_x, temp), hold on;
    plot(t_x, temp_fit)
    title('Ambient Temperature')
    xticks(t_x);
    xticklabels(time_labels);
end