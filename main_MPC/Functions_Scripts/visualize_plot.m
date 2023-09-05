

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
plot(squeeze(xx1(end,2,:))/par.E_bat_max)
plot(xSoC_N/par.E_bat_max);
legend('SoC(k)', 'SoC target DP(k)', 'SoC (predicted) (k+N)', 'SoC target (k+N)')
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

%%
figure

subplot(3,1,1)
plot(xx(1,:)*3.6,'k','LineWidth',1.), hold on
plot(linspace(mean(xx(1,:)*3.6),mean(xx(1,:)*3.6),length(xx(1,:))),'linewidth',1.)
plot(par.Route.max_v(1+iter_initial:iter_initial+size(xx,2))*3.6)
plot(par.v_DP*3.6,'g','LineWidth',1)
legend('velocity', 'mean velocity', 'max velocity','DP velocity')
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

subplot(3,1,2)
plot(xx(2,:)/par.E_bat_max,'k','LineWidth',1.0), hold on
plot(SoC_target(1+iter_initial:iter_initial+size(xx,2))/par.E_bat_max,'g','LineWidth',1.0)
plot(squeeze(xx1(end,2,:))/par.E_bat_max)
plot(xSoC_N/par.E_bat_max);
legend('SoC(k)', 'DP SoC target (k)', 'SoC (predicted) (k+N)', 'SoC target (k+N)')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[-]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,1,3)
plot(u_cl(:,1),'k','linewidth',.5), hold on
plot(linspace(mean(u_cl(:,1)),mean(u_cl(:,1)),length(u_cl(:,1))),'linewidth',1.)
plot(u_cl(:,2),'r','linewidth',.5)
plot(par.P_mot_el_DP,'g','LineWidth',1)
legend('input P_{el}', 'mean input P_{el}','input P_{br}','DP input P_{el}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[W]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

%% 
figure

plot(xSoC_diff/par.E_bat_max), hold on
plot(xS/par.E_bat_max)
legend('SoC(k+N)-SoC_{target}(k+N)', 'SoC_{deficit} (slack variable)')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[-]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*h+par.s_initial)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

%% 
