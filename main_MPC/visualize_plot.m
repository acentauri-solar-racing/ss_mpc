figure


subplot(3,3,1)
plot(xx(1,:)*3.6), hold on
plot(linspace(mean(xx(1,:)*3.6),mean(xx(1,:)*3.6),length(xx(1,:))),'linewidth',1.)
plot(par.Route.max_v(1+iter_initial:iter_initial+size(xx,2))*3.6)
legend('velocity', 'mean velocity', 'max velocity')
xlim([0 par.s_tot/par.s_step])
ylim([-10; 150])
xlabel('[km]')
ylabel('[km/h]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,4)
plot(xx(2,:)/par.E_bat_max), hold on
plot(SoC_target(1+iter_initial:iter_initial+size(xx,2))/par.E_bat_max)
legend('SoC', 'SoC target')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[-]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,7)
plot(u_cl(:,1),'k','linewidth',1.), hold on
plot(linspace(mean(u_cl(:,1)),mean(u_cl(:,1)),length(u_cl(:,1))),'linewidth',1.)
plot(u_cl(:,2),'r','linewidth',1.)
legend('input P_{el}', 'mean input P_{el}','input P_{br}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[W]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,2)
plot(simvar.alpha(1+iter_initial:iter_initial+size(xx,2)))
legend('\alpha')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[rad]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,5)
plot(xG)
legend('G')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[w/m^2]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,8)
plot(xtheta)
legend('\theta_{amb}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[°C]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);


subplot(3,3,3)
plot(xwind_front)
legend('v_{wind,front}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[m/s]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,6)
plot(xwind_side)
legend('v_{wind,side}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[m/s]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);



% figure
% 
% subplot(4,1,1)
% 
% plot(xx(3,2:end), u_cl(:,1)/par.P_el_max,'k','linewidth',1.5)
% legend('control input P_{el} (normalized)')
% xlabel('[s]')
% ylabel('[-]')
% title('time domain')
% 
% subplot(4,1,2)
% plot(xx(3,:), xx(1,:)*3.6), hold on
% plot(xx(3,:), par.Route.max_v(1:size(xx,2))*3.6)
% legend('velocity', 'max velocity')
% xlabel('[s]')
% ylabel('[km/h]')
% 
% 
% subplot(4,1,3)
% plot(xx(3,:), xx(2,:)/par.E_bat_max), hold on
% plot(xx(3,:), SoC_target(1:size(xx,2))/par.E_bat_max)
% legend('SoC', 'SoC target')
% xlabel('[s]')
% ylabel('[-]')
% 
% subplot(4,1,4)
% plot(xx(3,:), simvar.alpha(1:size(xx,2)))
% legend('\alpha')
% xlabel('[s]')
% ylabel('[rad]')
