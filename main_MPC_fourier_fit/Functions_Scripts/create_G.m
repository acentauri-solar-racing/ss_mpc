%% DP

% Define the length of the vector
vector_length = 540;

% Initialize the vector with a constant value of 800
G_MPC = 1000 * ones(1, vector_length);
G_real = 1000 * ones(1, vector_length);
G_DP = 1000 * ones(1, vector_length);

% Segment 1: Linearly increasing value from 300 to 400
G_MPC(1:21) = linspace(300, 400, 21);
G_real(1:21) = linspace(280, 350, 21);
G_DP(1:21) = linspace(300,400,21);

% Segment 2: Linearly increasing value from 400 to 600
G_MPC(21:91) = linspace(400, 600, 71);
G_real(21:81) = linspace(350, 570, 61);
G_DP(21:91) = linspace(400, 650, 71);

% Segment 3: Sinusoidal value between 550 and 650
segment_length = 30;
G_MPC(91:111) = 600; %+ 30 * sin(linspace(0, 4 * pi, 21));
G_real(82:111) = 570; %+ 30* sin(linspace(0, 8 * pi, segment_length));
G_DP(91:111) = 650;
% Segment 4: Constant value of 600

% Segment 5: Linearly increasing value to 800
G_MPC(112:211) = linspace(G_MPC(111), 500, 100);
G_real(112:181) = linspace(G_real(111), 490, 70);
G_DP(112:211) = linspace(G_DP(111), 475, 100);

G_MPC(211:241) = linspace(G_MPC(211), 1000, 31);
G_real(181:261) = linspace(G_real(181), 1000, 81);
G_DP(211:231) = linspace(G_DP(211), 1000, 21);

% Optional: Save the vector to a file
% save('G.mat', 'G');

save('..\OnlineData\G_MPC', "G_MPC")
save('..\OnlineData\G_real', "G_real")
save('..\OfflineData\G_DP', "G_DP")

% Plot the vector (optional)
plot(G_MPC), hold on;
plot(G_real);
plot(G_DP);
title('Generated 1D Vector');
xlim([0 250])
xlabel('time [min]');
ylabel('G [W/m^2]');
legend('G_{MPC}', 'G_{real}', 'G_{DP}')

% Optional: Save the vector to a file
% save('G.mat', 'G');