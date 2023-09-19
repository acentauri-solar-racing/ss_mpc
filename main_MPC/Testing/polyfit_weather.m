G_y = load('..\OnlineData\WeatherIrradiance').irradiance.Gtotal(3500:3620);
t_x = linspace(0,length(G_y),length(G_y))';

tic
p = polyfit(t_x,G_y,4);
toc
G_y_fit = polyval(p,t_x);

figure
plot(t_x,G_y), hold on;
plot(t_x,G_y_fit)


syms G(t)
G(t) = p(1)*t^2+p(2)*t+p(3)