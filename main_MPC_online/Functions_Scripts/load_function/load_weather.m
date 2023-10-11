%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 04.10.2023
%% Description 
% This function loads weatherdata from a certain folder
% INPUT : 
% "weatherfolderpath": string, folder path that contains the csv weather
% data

% OUTPUT : 
% "timeline (time)": vector that contains the time column resolution of the data
% "cumdist (space)": vector that contains the time cumulative distance resolution of the data
% "G_data": Global irradiance data (matrix dimension "time"x"space")
% "fW_data": Front Wind data (matrix dimension "time"x"space")
% "sW_data": Side Wind data (matrix dimension "time"x"space")
% "temp_data": Temperature data (matrix dimension "time"x"space")

%% Import data from text file


function [weather] = load_weather(weather, t_machine, weatherfolderpath)
    % create weather struct
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
    G_table = readtable(weatherfolderpath +"globalIrradiance.csv", opts);
    fW_table = readtable(weatherfolderpath + "frontWind.csv", opts);
    sW_table = readtable(weatherfolderpath + "sideWind.csv", opts);
    temp_table = readtable(weatherfolderpath + "temperature.csv", opts);
    rho_table = readtable(weatherfolderpath + "airDensity.csv", opts);
  
%   temp_table = readtable("C:\Users\loito\Desktop\alphacentauri_shared\ss_mpc\main_MPC_v3_polynomial_weather\OnlineData\20230926_090342_SF\preprocess\temperature.csv", opts);

%     G_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\globalIrradiance.csv", opts);
%     fW_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\frontWind.csv", opts);
%     sW_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\sideWind.csv", opts);
%     temp_table = readtable("G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\20231001_155728_SF\preprocess\temperature.csv", opts);
    
    timeline_string = cell2mat(table2array(G_table(2:end,1)));
    % get hour in "HHMM" format
    timeline_string = timeline_string(:,[12,13,15,16]);
    weather.timeline_string = timeline_string;
    weather.timeline = [];
    
    % get hour in seconds
    for i = 1:length(timeline_string)
        weather.timeline(i) = str2num(timeline_string(i,1))*10*60*60 + str2num(timeline_string(i,2))*60*60 + str2num(timeline_string(i,3))*10*60 + str2num(timeline_string(i,4))*60;
    end

    weather.timeline = weather.timeline';
    weather.cumdist = table2array(G_table(1,2:end));
    % Step 1: Create a logical index to find rows with 0
    zero_rows = any(weather.timeline == 0, 2);    
    % Step 2: Find the first row where 0 appears
    idx_zero = find(zero_rows, 1);
    
    % keep weather data until midnight
    if ~isempty(idx_zero)
        weather.G_data =  cellfun(@str2num,table2array(G_table(2:idx_zero, 2:end)));
        weather.fW_data = cellfun(@str2num,table2array(fW_table(2:idx_zero, 2:end)));
        weather.sW_data = cellfun(@str2num,table2array(sW_table(2:idx_zero, 2:end)));
        weather.temp_data = cellfun(@str2num,table2array(temp_table(2:idx_zero, 2:end)));
        weather.rho_data = cellfun(@str2num,table2array(rho_table(2:idx_zero, 2:end)));
        weather.timeline = weather.timeline(1:idx_zero-1);
    else  
        weather.G_data =  cellfun(@str2num,table2array(G_table(2:end, 2:end)));
        weather.fW_data = cellfun(@str2num,table2array(fW_table(2:end, 2:end)));
        weather.sW_data = cellfun(@str2num,table2array(sW_table(2:end, 2:end)));
        weather.temp_data = cellfun(@str2num,table2array(temp_table(2:end, 2:end)));
        weather.rho_data = cellfun(@str2num,table2array(rho_table(2:end, 2:end)));
    end
    
    % Get data; e.g. if you are at 9.35, get data from 9.30
    idx_overtime = find(weather.timeline > t_machine-15*60-1);
    weather.timeline = weather.timeline(idx_overtime);
    weather.G_data = weather.G_data(idx_overtime);
    weather.fW_data = weather.fW_data(idx_overtime);
    weather.sW_data = weather.sW_data(idx_overtime);
    weather.rho_data = weather.rho_data(idx_overtime);
    weather.temp_data = weather.temp_data(idx_overtime);


    %% Clear temporary variables
    clear opts
end