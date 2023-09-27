
% Example usage:
% Define the initial guess and the number of waves you want to use
numWavesToFit = 50; % Change this variable to adjust the number of waves

% Generate example data
% xData = linspace(0, 2 * pi, 100);
% yData = 2 * sin(2 * xData) + 0.5 * sin(3 * xData + pi/4) + 0.2 * sin(4 * xData - pi/3) + 0.1 * randn(size(xData));
t_start = 9515
time = 15
yData = load('..\OnlineData\WeatherIrradiance').irradiance.Gtotal(t_start:t_start+time)';
xData = linspace(0,length(yData),length(yData));

tic
% Initialize parameters for each wave in a cell array
initialGuess = cell(numWavesToFit, 1);

% Initialize parameters for all waves in a single vector
initialGuess = 1000*rand(1, 3 * numWavesToFit); % [A1, B1, C1, A2, B2, C2, ...]

% Create an options structure
options = optimoptions('lsqcurvefit', 'MaxFunctionEvaluations', 100000000);

% Use lsqcurvefit to fit the model to your data
paramsFit = lsqcurvefit(@(params, x) flexibleModel(params, x, numWavesToFit), initialGuess, xData, yData);

% Generate the fitted curve using the fitted parameters
yFit = flexibleModel(paramsFit, xData, numWavesToFit);
toc

% Plot the original data and the fitted curve
figure;
plot(xData, yData, 'o', xData, yFit, '-');
legend('Data', 'Fitted Curve');
xlabel('x');
ylabel('y');
title('Curve Fitting with Variable Number of Waves');

% Define the function to generate the model
function y = flexibleModel(params, x, numWaves)
    y = zeros(size(x));

    for i = 1:numWaves
        A = params((i - 1) * 3 + 1);
        B = params((i - 1) * 3 + 2);
        C = params((i - 1) * 3 + 3);
        y = y + A * sin(B * x + C);
    end
end