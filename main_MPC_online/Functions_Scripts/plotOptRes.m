function plotOptRes(OptRes, par)
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
    last100Data = max(1, numel(OptRes.xx(1, :)) - 50 + 1):numel(OptRes.xx(1, :));
    plot(last100Data, OptRes.xx(1, last100Data) * 3.6, 'LineWidth', 1.5, 'Color', 'b'), hold on;
    plot(last100Data, par.route.max_v(par.iter_initial + last100Data) * 3.6, 'LineWidth', 0.5, 'Color', 'r');
    title('velocity'); % Set a title
    ylim([0 140]);
    legend('predicted', 'max vel');
    xlabel('[km]'); % Label x-axis
    ylabel('[km/h]'); % Label y-axis

    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x * par.s_step + par.s_0) / 1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    % Update the second subplot (OptRes.xx(2, :))
    subplot(3, 1, 2);
    plot(last100Data, OptRes.xx(2, last100Data) / par.E_bat_max * 100, 'LineWidth', 1.5, 'Color', 'b'), hold on;
    plot(last100Data, par.E_bat_target_DP(par.iter_initial + last100Data) / par.E_bat_max * 100, 'LineWidth', 0.5, 'Color', 'r');
    ylim([OptRes.xx(2, end) / par.E_bat_max * 100 - 2.5 OptRes.xx(2, end) / par.E_bat_max * 100 + 2.5]);
    legend('predicted', 'target');
    title('SoC'); % Set a title
    xlabel('[km]'); % Label x-axis
    ylabel('[%]'); % Label y-axis

    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x * par.s_step + par.s_0) / 1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    % Update the third subplot (OptRes.u_cl(:,1))
    subplot(3, 1, 3);
    plot(last100Data, OptRes.u_cl(last100Data, 1), 'LineWidth', 1.5, 'Color', 'b'); % Plot OptRes.u_cl(:,1)
    title('Electric Motor Input'); % Set a title
    xlabel('[km]'); % Label x-axis
    ylabel('[W]'); % Label y-axis
    ylim([-5000 5000]);

    % Calculate the new tick positions and labels
    new_tick_positions = get(gca, 'XTick');
    new_tick_labels = arrayfun(@(x) num2str((x * par.s_step + par.s_0) / 1000), new_tick_positions, 'UniformOutput', false);
    % Set the new tick positions and labels
    set(gca, 'XTick', new_tick_positions);
    set(gca, 'XTickLabel', new_tick_labels);

    drawnow; % Refresh the figure window
end

