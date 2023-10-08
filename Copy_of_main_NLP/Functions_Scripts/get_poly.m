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
% G_1, G_2, G_3, G_4: second order polinomial fit parameters, G_1*t^3 + G_2*t^2 + G_3*t + G_4

%%
function [G_1, G_2, G_3, G_4] = get_poly(weather_data, t_0, t_f)
    
    G_y = weather_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1);
    t_x = linspace(floor(t_0/60/15),ceil(t_f/60/15),length(G_y))';

    p = polyfit(t_x,G_y,3);

    G_1 = p(1);
    G_2 = p(2);
    G_3 = p(3);
    G_4 = p(3);

end