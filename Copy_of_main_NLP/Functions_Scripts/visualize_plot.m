function visual = visualize_plot(par, OptResNLP)
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
    figure

end