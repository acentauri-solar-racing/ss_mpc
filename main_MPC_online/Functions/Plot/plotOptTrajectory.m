function plotOptTrajectory(OptRes, par, s, x0)
    % Check if it's the first iteration and create the figure and subplots
    if par.iter_mpc == 0
        figure; % Create a new figure

        % Create the first subplot for OptRes.xx(1, :)
        subplot(3, 1, 1);
        hold on;

        % Create the second subplot for OptRes.xx(2, :)
        subplot(3, 1, 2);
        hold on;

        % Create the third subplot for OptRes.u_cl(:,1)
        subplot(3, 1, 3);
        hold on;
    end

    % Clear the current subplot
    clf;

    % Update the first subplot (OptRes.xx(1, :))
    subplot(3, 1, 1);
    plot(1,x0(1)*3.6, 'Color','green', 'Marker','square', 'MarkerSize', 8), hold on;
    plot(OptRes.trajectory(2:end,1)*3.6, 'LineWidth', 1.5, 'Color', 'b')
    plot(par.route.max_v(par.iter_initial+par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)*3.6, 'LineWidth', 0.5, 'Color', 'r')
    title('velocity'); % Set a title
    ylim([40 140]);
    legend('car', 'predicted', 'max');
    xlabel('[km]'); % Label x-axis
    ylabel('[km/h]'); % Label y-axis

    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x * par.s_step + s) / 1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    % Update the second subplot (OptRes.xx(2, :))
    subplot(3, 1, 2);
    plot(1,x0(2) / par.E_bat_max * 100, 'Color','green', 'Marker','square', 'MarkerSize', 8), hold on;
    plot(OptRes.trajectory(2:end,2)/ par.E_bat_max * 100, 'LineWidth', 1.5, 'Color', 'b')
    plot(par.E_bat_target_DP(par.iter_initial +par.iter_mpc+1:par.iter_initial+par.N+1+par.iter_mpc)/ par.E_bat_max * 100, 'LineWidth', 0.5, 'Color', 'r')
    ylim([OptRes.xx(2, end) / par.E_bat_max * 100 - 2.5 OptRes.xx(2, end) / par.E_bat_max * 100 + 2.5]);
    legend('car', 'predicted', 'target');
    title('SoC'); % Set a title
    xlabel('[km]'); % Label x-axis
    ylabel('[%]'); % Label y-axis

    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x * par.s_step + s) / 1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    % Update the third subplot (OptRes.u_cl(:,1))
    subplot(3, 1, 3);
    plot(OptRes.u_trajectory(:,1), 'LineWidth', 1.5, 'Color', 'b')
    title('Electric Motor Input'); % Set a title
    xlabel('[km]'); % Label x-axis
    ylabel('[W]'); % Label y-axis
    ylim([-2500 2500]);
    mean_input = mean(OptRes.u_trajectory(:,1));
    yline(mean_input,'-', ['mean = ', num2str(mean_input)] ,'Color','g');

    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x * par.s_step + s) / 1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);


    drawnow; % Refresh the figure window
end

