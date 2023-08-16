figure


subplot(3,2,1)
stairs(xx(1,:)*3.6), hold on
stairs(par.Route.max_v(1+iter_initial:iter_initial+size(xx,2))*3.6)
legend('velocity', 'max velocity')
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

subplot(3,2,3)
stairs(xx(2,:)/par.E_bat_max), hold on
stairs(SoC_target(1+iter_initial:iter_initial+size(xx,2))/par.E_bat_max)
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

subplot(3,2,5)
stairs(u_cl(:,1),'k','linewidth',1.), hold on
stairs(u_cl(:,2),'r','linewidth',1.)
legend('input P_{el}', 'input P_{br}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[W]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,2,2)
stairs(simvar.alpha(1+iter_initial:iter_initial+size(xx,2)))
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

subplot(3,2,4)
stairs(simvar.G(1+iter_initial:iter_initial+size(xx,2)))
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

% figure
% 
% subplot(4,1,1)
% 
% stairs(xx(3,2:end), u_cl(:,1)/par.P_el_max,'k','linewidth',1.5)
% legend('control input P_{el} (normalized)')
% xlabel('[s]')
% ylabel('[-]')
% title('time domain')
% 
% subplot(4,1,2)
% stairs(xx(3,:), xx(1,:)*3.6), hold on
% stairs(xx(3,:), par.Route.max_v(1:size(xx,2))*3.6)
% legend('velocity', 'max velocity')
% xlabel('[s]')
% ylabel('[km/h]')
% 
% 
% subplot(4,1,3)
% stairs(xx(3,:), xx(2,:)/par.E_bat_max), hold on
% stairs(xx(3,:), SoC_target(1:size(xx,2))/par.E_bat_max)
% legend('SoC', 'SoC target')
% xlabel('[s]')
% ylabel('[-]')
% 
% subplot(4,1,4)
% stairs(xx(3,:), simvar.alpha(1:size(xx,2)))
% legend('\alpha')
% xlabel('[s]')
% ylabel('[rad]')
