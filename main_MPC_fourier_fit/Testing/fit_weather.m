%G_y = load('..\OnlineData\G_MPC').G_MPC(75:100);
y = load('..\OnlineData\WeatherIrradiance').irradiance.Gtotal(3530+225:3530+235 )';
x = linspace(0,length(y),length(y));
% Initial parameter guesses (adjust these as needed)
numWaves = 3;  % Number of waves in the superposition model

% Define your custom model function for a variable number of waves
model = @(params, x) sum(params(1:numWaves) .* sin(params(numWaves+1:2*numWaves) * x + params(2*numWaves+1:end)));

% Initialize initial guesses for the parameters
initialGuesses = zeros(1, 3 * numWaves);

% Fit the model to your data using lsqcurvefit
options = optimoptions('lsqcurvefit', 'Display', 'off');
params = lsqcurvefit(model, initialGuesses, x, y);

% Extract the coefficients of the fitted model
amplitudes = params(1:numWaves);
frequencies = params(numWaves+1:2*numWaves);
phases = params(2*numWaves+1:end);

% Plot the data and the fitted curve
fittedCurve = zeros(size(x));
for i = 1:numWaves
    fittedCurve = fittedCurve + amplitudes(i) * sin(frequencies(i) * x + phases(i));
end

plot(x, y, 'o', x, fittedCurve, '-');

% Display the coefficients of the fitted model
disp('Coefficients:');
for i = 1:numWaves
    disp(['Amplitude (a', num2str(i), '): ', num2str(amplitudes(i))]);
    disp(['Frequency (b', num2str(i), '): ', num2str(frequencies(i))]);
    disp(['Phase (c', num2str(i), '): ', num2str(phases(i))]);
end
%%
% t_x = linspace(0,length(G_y),length(G_y))';
% 
% tic
% p = polyfit(t_x,G_y,2);
% toc
% G_y_fit = polyval(p,t_x);
% 
% figure
% plot(t_x,G_y), hold on;
% plot(t_x,G_y_fit)
% legend('real','fit')
% 
% 
% syms G(t)
% G(t) = p(1)*t^2+p(2)*t+p(3)

% %%
% G_y = load('..\OnlineData\WeatherIrradiance').irradiance.Gtotal(3500:3600);
% t_min = linspace(0,length(G_y),length(G_y))';
% t_sec = linspace(0,length(G_y)*60,length(G_y)*60+1)';
% 
% G_y_sec = interp1(t_min,G_y,t_sec);
% 
% G_y_min = G_y_sec(1:60:length(G_y_sec)-1);
% 
% p = polyfit(t_min,G_y_min,4);
% 
% 
% plot(t_sec,G_y_sec)