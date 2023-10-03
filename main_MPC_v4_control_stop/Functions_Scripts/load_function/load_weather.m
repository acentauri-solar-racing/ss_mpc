
%% Import data from text file


function [time, cumdist, G_data, fW_data, sW_data, temp_data] = load_weather()
    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 11);
    
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    % opts.VariableNames = ["time", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11"];
    % opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    
    % % Specify variable properties
    % opts = setvaropts(opts, "time", "WhitespaceRule", "preserve");
    % opts = setvaropts(opts, "time", "EmptyFieldRule", "auto");
    
    % Import the data
    G_table = readtable("C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_polynomial_weather\OnlineData\20230926_090342_SF\preprocess\globalIrradiance.csv", opts);
    fW_table = readtable("C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_polynomial_weather\OnlineData\20230926_090342_SF\preprocess\frontWind.csv", opts);
    sW_table = readtable("C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_polynomial_weather\OnlineData\20230926_090342_SF\preprocess\sideWind.csv", opts);
    temp_table = readtable("C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_polynomial_weather\OnlineData\20230926_090342_SF\preprocess\temperature.csv", opts);
%         G_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\globalIrradiance.csv", opts);
%     fW_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\frontWind.csv", opts);
%     sW_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\sideWind.csv", opts);
%     temp_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\temperature.csv", opts);
    time = cell2mat(table2array(G_table(2:end,1)));
    cumdist = table2array(G_table(1,2:end));
    
    G_data =  cellfun(@str2num,table2array(G_table(2:end, 2:end)));
    fW_data = cellfun(@str2num,table2array(fW_table(2:end, 2:end)));
    sW_data = cellfun(@str2num,table2array(sW_table(2:end, 2:end)));
    temp_data = cellfun(@str2num,table2array(temp_table(2:end, 2:end)));
    
    
    %% Clear temporary variables
    clear opts
end