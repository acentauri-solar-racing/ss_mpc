t_0 = 0;
t_f = 60*15*4*4;

G_y =   temp_mpc(1+floor(t_0/60/15):1+ceil(t_f/60/15));
t_x = linspace(floor(t_0/60/15),ceil(t_f/60/15),length(G_y))';
tic
p = polyfit(t_x,G_y,2);
toc
G_y_fit = polyval(p,t_x);


%%
figure
plot(t_x,G_y), hold on;
plot(t_x,G_y_fit)


syms G(t)
G(t) = p(1)*t^2+p(2)*t+p(3)

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
%%
