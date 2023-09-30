figure

subplot(3,3,1)
plot(OptRes.xx(1,:)*3.6), hold on
plot(linspace(mean(OptRes.xx(1,:)*3.6),mean(OptRes.xx(1,:)*3.6),length(OptRes.xx(1,:))),'linewidth',1.)
plot(par.Route.max_v(1+par.iter_initial:par.iter_initial+size(OptRes.xx,2))*3.6)

title('velocity')
legend('velocity', 'mean velocity', 'max velocity')
xlim([0 par.s_tot/par.s_step])
ylim([-10; 150])
xlabel('[km]')
ylabel('[km/h]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,4)
plot(OptRes.xx(2,:)/par.E_bat_max), hold on
plot(par.E_bat_target_DP(1+par.iter_initial:par.iter_initial+size(OptRes.xx,2))/par.E_bat_max)
plot(squeeze(OptRes.xx1(end,2,:))/par.E_bat_max)
plot(OptRes.xSoC_N/par.E_bat_max);
title('State of Charge')
legend('SoC(k)', 'SoC target DP(k)', 'SoC (predicted) (k+N)', 'SoC target (k+N)')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[-]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,7)
plot(OptRes.u_cl(:,1),'k','linewidth',1.), hold on
plot(linspace(mean(OptRes.u_cl(:,1)),mean(OptRes.u_cl(:,1)),length(OptRes.u_cl(:,1))),'linewidth',1.)
plot(OptRes.u_cl(:,2),'r','linewidth',1.)
title('Control Input')
legend('input P_{el}', 'mean input P_{el}','input P_{br}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[W]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,2)
plot(par.Route.incl(1+par.iter_initial:par.iter_initial+size(OptRes.xx,2)))
title('Inclination')
legend('\alpha')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[rad]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,5)
plot(OptRes.xG)
title('Solar Irradiation')
legend('G')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[w/m^2]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,8)
plot(OptRes.xtheta)
title('Amb. Temperature')
legend('\theta_{amb}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[°C]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);


subplot(3,3,3)
plot(OptRes.xwind_front)
title('vel. front wind')
legend('v_{wind,front}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[m/s]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,6)
plot(OptRes.xwind_side)
title('vel. side wind')
legend('v_{wind,side}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[m/s]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,3,9)
plot(OptRes.xx(3,:))
title('time')
legend('time')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[s]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

%%
figure

subplot(3,1,1)
plot(OptRes.xx(1,:)*3.6,'k','LineWidth',1.), hold on
plot(linspace(mean(OptRes.xx(1,:)*3.6),mean(OptRes.xx(1,:)*3.6),length(OptRes.xx(1,:))),'linewidth',1.)
plot(par.Route.max_v(1+par.iter_initial:par.iter_initial+size(OptRes.xx,2))*3.6)
plot(par.v_DP*3.6,'g','LineWidth',1)
title('velocity')
legend('velocity', 'mean velocity', 'max velocity','DP velocity')
xlim([0 par.s_tot/par.s_step])
ylim([-10; 150])
xlabel('[km]')
ylabel('[km/h]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,1,2)
plot(OptRes.xx(2,:)/par.E_bat_max,'k','LineWidth',1.0), hold on
plot(par.E_bat_target_DP(1+par.iter_initial:par.iter_initial+size(OptRes.xx,2))/par.E_bat_max,'g','LineWidth',1.0)
plot(squeeze(OptRes.xx1(end,2,:))/par.E_bat_max)
plot(OptRes.xSoC_N/par.E_bat_max);
title('State of Charge')
legend('SoC(k)', 'DP SoC target (k)', 'SoC (predicted) (k+N)', 'SoC target (k+N)')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[-]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

subplot(3,1,3)
plot(OptRes.u_cl(:,1),'k','linewidth',.5), hold on
plot(linspace(mean(OptRes.u_cl(:,1)),mean(OptRes.u_cl(:,1)),length(OptRes.u_cl(:,1))),'linewidth',1.)
plot(OptRes.u_cl(:,2),'r','linewidth',.5)
plot(par.P_mot_el_DP,'g','LineWidth',1)
title('Control Input')
legend('input P_{el}', 'mean input P_{el}','input P_{br}','DP input P_{el}')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[W]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

%% 
figure

plot(OptRes.xSoC_diff/par.E_bat_max), hold on
plot(OptRes.xS/par.E_bat_max)
legend('SoC(k+N)-SoC_{target}(k+N)', 'SoC_{deficit} (slack variable)')
xlim([0 par.s_tot/par.s_step])
xlabel('[km]')
ylabel('[-]')

% Calculate the new tick positions and labels
new_tick_positions = get(gca, 'XTick');
new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
% Set the new tick positions and labels
set(gca, 'XTick', new_tick_positions);
set(gca, 'XTickLabel', new_tick_labels);

%% 
