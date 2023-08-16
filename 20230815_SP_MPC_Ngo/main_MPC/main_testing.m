clear all
close all
clc

%%

% Assuming you have two datasets: data1 and data2
xx100 = load('100.mat').xx;
xx100euler = load('100euler.mat').xx;
xx50 = load('50.mat').xx;

common_length = 200; % Choose a common length for resampling

% Create the x-axis for the first dataset, size(xx100,2)
x1 = 0:100:(100 * (size(xx100,2) - 1));

% Create the x-axis for the second dataset
x2 = 0:100:(100 * (size(xx100euler,2) - 1));

figure;
plot(x1, xx100(1,:), 'b', x2, xx100euler(1,:), 'r');
legend('Dataset 1', 'Dataset 2');
xlabel('Distance (m)');
ylabel('Data Value');
title('Comparison of Datasets');
