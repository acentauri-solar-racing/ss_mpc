function visual = visualize_plot_BWSC(par, weather, OptResNLP)
%% NLP only States/Input
    figure
    subplot(4,2,1)
    plot(OptResNLP.xx(:,1)*3.6, 'LineStyle','-','linewidth',1.5), hold on
    plot(par.route.max_v(1+par.iter_initial:par.iter_initial+par.N+1)*3.6,'-.', 'LineWidth', 0.5, 'Color', 'k')
    
    title('velocity')
    legend('nlp', 'max velocity')
    xlim([0 par.s_tot/par.s_step])
    ylim([0; 150])
    xlabel('[km]')
    ylabel('[km/h]')
    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);
    
    subplot(4,2,3)
    plot(OptResNLP.xx(:,2)/par.E_bat_max, 'LineStyle','-','linewidth',1.5), hold on
    plot(par.E_bat_target_DP(1+par.iter_initial:par.iter_initial+par.N+1)/par.E_bat_max)
    
    title('SoC')
    legend('nlp', 'target')
    xlim([0 par.s_tot/par.s_step])
    xlabel('[km]')
    ylabel('[-]')
    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    subplot(4,2,5)
    plot(OptResNLP.u_cl(:,1),'linewidth',1.5, 'LineStyle','-'), hold on
    plot(OptResNLP.u_cl(:,2),'linewidth',0.1,'color','r', 'LineStyle','--')
    
    title('Control Input')
    legend('nlp', 'nlp (brake)')
    xlim([0 par.s_tot/par.s_step])
    xlabel('[km]')
    ylabel('[W]')
    
    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    %% Weather
    t = OptResNLP.xx(:,3);
    t_15min = OptResNLP.xx(:,3)/60/15; 
    t_0 = par.t_0;
    t_f = OptResNLP.xx(end,3);
    t_15min = linspace(floor(t_0/60/15),ceil(t_f/60/15),length(weather.G_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1)))';

    G= par.G_1*(t/60/15).^3 + par.G_2*(t/60/15).^2 + par.G_3*(t/60/15) + par.G_4;
    fW= par.fW_1*(t/60/15).^3 + par.fW_2*(t/60/15).^2 + par.fW_3*(t/60/15) + par.fW_4;
    rho= par.rho_1*(t/60/15).^3 + par.rho_2*(t/60/15).^2 + par.rho_3*(t/60/15) + par.rho_4;
    temp= par.temp_1*(t/60/15).^3 + par.temp_2*(t/60/15).^2 + par.temp_3*(t/60/15) + par.temp_4;

    G_real = interp1(t_15min, weather.G_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx(:,3)/60/15);
    fW_real = interp1(t_15min, weather.fW_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx(:,3)/60/15);
    rho_real = interp1(t_15min, weather.rho_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx(:,3)/60/15);
    temp_real = interp1(t_15min, weather.temp_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx(:,3)/60/15);



    
    subplot(4,2,2)
    plot(G), hold on;
    plot(G_real)
    title('Solar Irradiation')
    legend('fit', 'real')
    xlim([0 par.s_tot/par.s_step])
    xlabel('[km]')
    ylabel('[W/m^2]')
    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    subplot(4,2,4)
    plot(fW), hold on;
    plot(fW_real)
    title('Frontal Wind')
    legend('fit', 'real')
    xlim([0 par.s_tot/par.s_step])
    xlabel('[km]')
    ylabel('[m/s]')
    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    subplot(4,2,6)
    plot(rho), hold on;
    plot(rho_real)
    title('Air Density')
    legend('fit', 'real')
    xlim([0 par.s_tot/par.s_step])
    xlabel('[km]')
    ylabel('[kg/m^3]')
    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    subplot(4,2,8)
    plot(temp), hold on;
    plot(temp_real)
    title('Temperature')
    legend('fit', 'real')
    xlim([0 par.s_tot/par.s_step])
    xlabel('[km]')
    ylabel('[°C]')
    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x*par.s_step+par.s_0)/1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);
end