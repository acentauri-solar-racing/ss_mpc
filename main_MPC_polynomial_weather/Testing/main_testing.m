clear all
close all
clc

%%
% 
% % Assuming you have two datasets: data1 and data2
% xx100 = load('100.mat').xx;
% xx100euler = load('100euler.mat').xx;
% xx50 = load('50.mat').xx;
% 
% common_length = 200; % Choose a common length for resampling
% 
% % Create the x-axis for the first dataset, size(xx100,2)
% x1 = 0:100:(100 * (size(xx100,2) - 1));
% 
% % Create the x-axis for the second dataset
% x2 = 0:100:(100 * (size(xx100euler,2) - 1));
% 
% figure;
% plot(x1, xx100(1,:), 'b', x2, xx100euler(1,:), 'r');
% legend('Dataset 1', 'Dataset 2');
% xlabel('Distance (m)');
% ylabel('Data Value');
% title('Comparison of Datasets');

%% 
% 271 - 1231
% load WeatherIrradiance.mat
% G = irradiance.Gtotal(271:1230);
% t = linspace(60*4.5,60*20.5,60*(20.5-4.5));
% G_int = griddedInterpolant(t,G);
% 
% % min = seconds/60
% vG = G_int(40235.43/60)

% theta_2d = load("WeatherData.mat").temperature.tempMean(:,5:14); % from 8.30 to 17.30
% theta_dist = load("WeatherData.mat").temperature.dist;
% theta_time = linspace(0,9,10);
% 
% [theta_x, theta_t] = ndgrid(theta_dist,theta_time);
% theta_int = griddedInterpolant(theta_x,theta_t,theta_2d);
% 
% theta_res = theta_int({linspace(0,100000,5),linspace(0,5,6)})
step = 10;
step_DP = 10000;
E_bat_target_DP_raw = load('Full_Race_20230607_9h_30min.mat').OptRes.states.E_bat*3600; % [Wh = W * 3600s = J]
E_bat_target_DP = interp1(linspace(0,300,300+1), E_bat_target_DP_raw, linspace(0,300,(300+1)*round(step_DP/step)));
