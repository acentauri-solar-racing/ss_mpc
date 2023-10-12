%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 10.11.2023
%% Description 
% This function save simulation data

% INPUT: 
% "OptRes" (struct): Simulation data stored struct
% "par" (struct): parameters struct
% "weather" (struct): weather info as struct
% "filename_date" (string): the date inserted as a string
% OUTPUT : 
% save files in the path 'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\MPC_'+filename_date+'\MPC_'+filename_date

function save_file = save_file(OptRes, par, weather, filename_date)
    % save mat file with OptRes, par and weather
    mkdir('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\MPC_'+filename_date)
    save('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\MPC_'+filename_date+'\MPC_'+filename_date, "OptRes","par","weather")
    writematrix([vertcat(OptRes.xdist, OptRes.xx)],'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\MPC_'+filename_date+'\MPC_'+filename_date+'_states.csv')
    writematrix([vertcat(OptRes.xdist(2:end), OptRes.u_cl')],'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\MPC_'+filename_date+'\MPC_'+filename_date+'_inputs.csv')
end
