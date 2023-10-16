%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 10.11.2023
%% Description 
% This function save simulation data

% INPUT: 
% "OptRes" (struct): Simulation data stored struct
% "par" (struct): parameters struct
% "weather" (struct): weather info as struct
% OUTPUT : 
% save files in the path 'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\MPC_'+par.filename+'\MPC_'+par.filename

function save_file = save_file_BWSC(OptRes, par, weather)

% Convert seconds to 'HH:MM:SS' format
timeStrings = datestr(seconds(OptRes.xx(:,3)), 'HH:MM:SS');
newTimeStrings = [];
for i = 1:length(timeStrings(:,1))
    newTimeStrings = [newTimeStrings;weather.ddmmyy_string,' ',timeStrings(i,:)];
end
state = table(newTimeStrings, OptRes.xdist,OptRes.xx(:,1)*3.6,OptRes.xx(:,2)/par.E_bat_max, 'VariableNames',{'time';'cumDistance';'velocity';'soc'});
input = table(newTimeStrings(2:end,:), OptRes.xdist(2:end),OptRes.u_cl(:,1),OptRes.u_cl(:,2), 'VariableNames',{'time';'cumDistance';'inputElectricMotorPower';'inputBrakePower'});

mkdir('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\'+par.filename+'_NLP')
save('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\'+par.filename+'_NLP\'+par.filename, "OptRes","par","weather")

writetable(state,'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\'+par.filename+'_NLP\'+par.filename+'_state_NLP.csv')
writetable(input,'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\'+par.filename+'_NLP\'+par.filename+'_input_NLP.csv')

%writematrix([vertcat(OptRes.xdist(2:end), OptRes.u_cl')],'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\'+par.filename+'_NLP\'+par.filename+'_inputs_NLP.csv')
end

%     % save mat file with OptRes, par and weather
%     currentDateTime = datetime('now');
%     timestamp = datestr(currentDateTime, 'yyyymmdd_HHMMSS');
%     directory = 'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\';
%     filename = [directory, timestamp, '_MPC'];
%     mkdir(filename);
%     filename = [filename, timestamp, '_MPC'];
%     save('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\BWSC_'+par.filename+'\BWSC_'+par.filename, "OptRes","par","weather")
%     writematrix([horzcat(OptRes.xdist, OptRes.xx)],'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\BWSC_'+par.filename+'\BWSC_'+par.filename+'_states.csv')
%     writematrix([horzcat(OptRes.xdist(2:end), OptRes.u_cl)],'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\BWSC_'+par.filename+'\BWSC_'+par.filename+'_inputs.csv')

