function visual = visualize_plot(par, weather, OptResNLP)
%% NLP only
    figure
    subplot(3,1,1)
    plot(OptResNLP.xx1(:,1)*3.6, 'LineStyle','-','linewidth',1.5), hold on
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
    
    subplot(3,1,2)
    
    plot(OptResNLP.xx1(:,2)/par.E_bat_max, 'LineStyle','-','linewidth',1.5), hold on
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
    subplot(3,1,3)
    
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

    %% 
    t = OptResNLP.xx1(:,3);
    t_15min = OptResNLP.xx1(:,3)/60/15; 
    t_0 = par.t_0;
    t_f = OptResNLP.xx1(end,3);
    t_15min = linspace(floor(t_0/60/15),ceil(t_f/60/15),length(weather.G_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1)))';

    G= par.G_1*(t/60/15).^3 + par.G_2*(t/60/15).^2 + par.G_3*(t/60/15) + par.G_4;
    fW= par.fW_1*(t/60/15).^3 + par.fW_2*(t/60/15).^2 + par.fW_3*(t/60/15) + par.fW_4;
    sW= par.sW_1*(t/60/15).^3 + par.sW_2*(t/60/15).^2 + par.sW_3*(t/60/15) + par.sW_4;
    temp= par.temp_1*(t/60/15).^3 + par.temp_2*(t/60/15).^2 + par.temp_3*(t/60/15) + par.temp_4;

    G_real = interp1(t_15min, weather.G_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx1(:,3)/60/15);
    fW_real = interp1(t_15min, weather.fW_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx1(:,3)/60/15);
    sW_real = interp1(t_15min, weather.sW_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx1(:,3)/60/15);
    temp_real = interp1(t_15min, weather.temp_data(1+floor(t_0/60/15):1+ceil(t_f/60/15),1),  OptResNLP.xx1(:,3)/60/15);



    figure
    subplot(4,1,1)
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

    subplot(4,1,2)
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

    subplot(4,1,3)
    plot(sW), hold on;
    plot(sW_real)
    title('Side Wind')
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

    subplot(4,1,4)
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