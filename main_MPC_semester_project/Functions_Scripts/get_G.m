%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 08.09.2023
%% Description 
% This function returns the polinomial fit parameters for the solar
% irradiation G, assuming it's only depending on the time

% INPUT: 
% "par": parameters struct
% "t_0": initial time in seconds of polinomial fit
% "t_f": final time in seconds of polinomial fit


% OUTPUT : 
% G_1, G_2, G_3: second order polinomial fit parameters, G_1*t^2 + G_2*t + G_3

%%
function [G_1, G_2, G_3] = get_G(G_data, t_0, t_f)
    G_y = G_data(1+round(t_0/60):1+round(t_f/60));
    t_x = linspace(round(t_0/60),round(t_0/60)+length(G_y),length(G_y))';

    p = polyfit(t_x,G_y,2);

    G_1 = p(1);
    G_2 = p(2);
    G_3 = p(3);
end