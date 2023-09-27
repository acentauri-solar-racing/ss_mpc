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
function G_params = get_G(n_waves, G_data, t_0, t_f)

    G_y = G_data(1+round(t_0/60):1+round(t_f/60))';
    t_x = linspace(round(t_0/60),round(t_0/60)+length(G_y),length(G_y))';
    
    
    % Initialize parameters for all waves in a single vector
    init_guess = 1000*rand(1, 3 * n_waves); % [A1, B1, C1, A2, B2, C2, ...]
    
    % Create an options structure
    options = optimoptions('lsqcurvefit', 'MaxFunctionEvaluations', 100000000, 'Algorithm','levenberg-marquardt');
    lb = [];
    ub = [];
    
    % Use lsqcurvefit to fit the model to your data
    G_params = lsqcurvefit(@(params, x) G_model(params, x, n_waves), init_guess, t_x, G_y,lb,ub,options);
    G_params = G_params';

end

function y = G_model(params, x, n_waves)
    y = zeros(size(x));

    for i = 1:n_waves
        A = params((i - 1) * 3 + 1);
        B = params((i - 1) * 3 + 2);
        C = params((i - 1) * 3 + 3);
        y = y + A * sin(B * x + C);
    end
end