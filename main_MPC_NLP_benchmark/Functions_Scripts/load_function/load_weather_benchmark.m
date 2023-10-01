function [t, G_nlp, G_mpc,  fW_nlp, fW_mpc, sW_nlp, sW_mpc, temp_nlp, temp_mpc] = load_weather_benchmark()


    %% irradiance
    t = [0:1:16]; % 1 = 15 min : 1 = 15  min : 16 = 240 min = 4 h
    G_nlp = -2.1*t.^2+62*t+420;
    G_nlp = G_nlp';
    G_mpc = -2*t.^2+55*t+470;
    G_mpc = G_mpc';
    
%     figure
%     plot(t,G_nlp), hold on
%     plot(t,G_mpc)
%     legend('nlp (reality)','mpc')
%     
%     % Calculate the new tick positions and labels
%     new_tick_positions = get(gca, 'XTick');
%     new_tick_labels = arrayfun(@(x) num2str(x*0.25), new_tick_positions, 'UniformOutput', false);
%     % Set the new tick positions and labels
%     set(gca, 'XTick', new_tick_positions);
%     set(gca, 'XTickLabel', new_tick_labels);
%     xlabel('time [h]')
%     ylabel('irradiance [W/m^2]')
    
    %% front wind
    
    fW_nlp = -0.02*t.^2+0.04*t+2.2;  
    fW_nlp = -0.00*t.^2+0.00*t+0.0;

    fW_nlp = fW_nlp';
    fW_mpc = -0*t.^0+0*t+0;
    fW_mpc = fW_mpc';

%     figure
%     plot(t,fW_nlp), hold on
%     plot(t,fW_mpc)
%     legend('nlp (reality)','mpc')
%     
%     % Calculate the new tick positions and labels
%     new_tick_positions = get(gca, 'XTick');
%     new_tick_labels = arrayfun(@(x) num2str(x*0.25), new_tick_positions, 'UniformOutput', false);
%     % Set the new tick positions and labels
%     set(gca, 'XTick', new_tick_positions);
%     set(gca, 'XTickLabel', new_tick_labels);
%     xlabel('time [h]')
%     ylabel('front wind velocity [m/s]')
    
    %% side wind
    
    sW_nlp = 0*t.^0+0*t+0;
    sW_nlp = sW_nlp';
    sW_mpc = 0*t.^0+0*t+0;
    sW_mpc = sW_mpc';

    %% temperature
    
    temp_nlp = 0.04*t.^2-0.06*t+25;
    temp_nlp = 0.0*t.^2-0.0*t+25;
    temp_nlp = temp_nlp';
    temp_mpc = 0.02*t.^2-0.04*t+27;
    temp_mpc = temp_mpc';
    
%     figure
%     plot(t,temp_nlp), hold on
%     plot(t,temp_mpc)
%     legend('nlp (reality)','mpc')
%     
%     % Calculate the new tick positions and labels
%     new_tick_positions = get(gca, 'XTick');
%     new_tick_labels = arrayfun(@(x) num2str(x*0.25), new_tick_positions, 'UniformOutput', false);
%     % Set the new tick positions and labels
%     set(gca, 'XTick', new_tick_positions);
%     set(gca, 'XTickLabel', new_tick_labels);
%     xlabel('time [h]')
%     ylabel('temperature [°C]')

end