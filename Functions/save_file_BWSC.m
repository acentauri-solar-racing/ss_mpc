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
    % save mat file with OptRes, par and weather
    mkdir('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\NLP_'+par.filename)
    save('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\NLP_'+par.filename+'\NLP_'+par.filename, "OptRes","par","weather")
    writematrix([horzcat(OptRes.xdist, OptRes.xx)],'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\NLP_'+par.filename+'\NLP_'+par.filename+'_states.csv')
    writematrix([horzcat(OptRes.xdist(2:end), OptRes.u_cl)],'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\MPC_optimal\NLP_'+par.filename+'\NLP_'+par.filename+'_inputs.csv')
end
