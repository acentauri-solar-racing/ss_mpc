%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 11.10.2023
%% Description 
% This function loads weatherdata from a certain folder; 
% IT'S IMPORTANT THAT "TIME MACHINE" > "WEATHER FORECAST START TIME", of
% the same day

% INPUT : 
% "t_0" [s]: actual time machine (time on the computer clock)
% "weather": weather struct
% "n_day": the day with respect the first day of data. n_day = 0 means you
% use the date of the first day until midnight the data; n_day = 1 would
% mean ignore the data of the first day and use the data of the day
% tomorrow, and so on

% OUTPUT : 
% "timeline (time)": vector that contains the time column resolution of the data
% "cumdist (space)": vector that contains the time cumulative distance resolution of the data
% "G_data": Global irradiance data (matrix dimension "time"x"space")
% "fW_data": Front Wind data (matrix dimension "time"x"space")
% "sW_data": Side Wind data (matrix dimension "time"x"space")
% "temp_data": Temperature data (matrix dimension "time"x"space")

%% Import data from text file


function [weather] = load_weather(weather, t_0, n_day)
    % get latest weather folder
    weatherfolderpath = get_latest_weather('G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\Forecast\');

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
    
    % Import the data as matrices
    G_table = readtable(weatherfolderpath +"globalIrradiance.csv", opts);
    fW_table = readtable(weatherfolderpath + "frontWind.csv", opts);
    sW_table = readtable(weatherfolderpath + "sideWind.csv", opts);
    temp_table = readtable(weatherfolderpath + "temperature.csv", opts);
    rho_table = readtable(weatherfolderpath + "airDensity.csv", opts);
    
    % Get time vector as string
    weather.fulldate_string = cell2mat(table2array(G_table(2:end,1)));
    % get hour in "HHMM" format
    timeline_string = weather.fulldate_string(:,[12,13,15,16]);
    weather.timeline_string = timeline_string;
    weather.timeline = [];
    
    % convert time vector string into seconds
    for i = 1:length(timeline_string)
        weather.timeline(i) = str2num(timeline_string(i,1))*10*60*60 + str2num(timeline_string(i,2))*60*60 + str2num(timeline_string(i,3))*10*60 + str2num(timeline_string(i,4))*60;
    end
    weather.timeline = weather.timeline';

    % get cumulative distance vector
    weather.cumdist_cell = table2array(G_table(1,2:end));
    % Preallocate the cumdist vector to match the size of cumdist_cell
    weather.cumdist = zeros(1, numel(weather.cumdist_cell));
    
    % Extract elements from the cell array, convert them to the appropriate data type (e.g., double), and store them in the vector
    for i = 1:numel(weather.cumdist_cell)
        weather.cumdist(1,i) = str2double(weather.cumdist_cell{1,i});
    end
    
    % Step 1: Create a logical index to find rows with 0
    zero_rows = any(weather.timeline == 0, 2);    
    % Step 2: Find the first row where 0 appears
    
    if n_day == 0;
        % keep weather data from start of datatime until midnight
        if ~isempty(zero_rows)
            idx_midnight = find(zero_rows, 1)-1;

            weather.G_data =  cellfun(@str2num,table2array(G_table(2:idx_midnight, 2:end)));
            weather.fW_data = cellfun(@str2num,table2array(fW_table(2:idx_midnight, 2:end)));
            weather.sW_data = cellfun(@str2num,table2array(sW_table(2:idx_midnight, 2:end)));
            weather.temp_data = cellfun(@str2num,table2array(temp_table(2:idx_midnight, 2:end)));
            weather.rho_data = cellfun(@str2num,table2array(rho_table(2:idx_midnight, 2:end)));
            weather.timeline = weather.timeline(1:idx_midnight-1);
            weather.timeline_string = weather.timeline_string(1:idx_midnight-1,:);
            weather.fulldate_string = weather.fulldate_string(1:idx_midnight-1,:);

        else  
            weather.G_data =  cellfun(@str2num,table2array(G_table(2:end, 2:end)));
            weather.fW_data = cellfun(@str2num,table2array(fW_table(2:end, 2:end)));
            weather.sW_data = cellfun(@str2num,table2array(sW_table(2:end, 2:end)));
            weather.temp_data = cellfun(@str2num,table2array(temp_table(2:end, 2:end)));
            weather.rho_data = cellfun(@str2num,table2array(rho_table(2:end, 2:end)));
            weather.timeline = weather.timeline(1:idx_midnight-1);
            weather.timeline_string = weather.timeline_string(1:idx_midnight-1,:);
            weather.fulldate_string = weather.fulldate_string(1:idx_midnight-1,:);
        end
    else
        idx_midnight = find(zero_rows, n_day+1);
        idx_start = idx_midnight(end-1);
        idx_end = idx_midnight(end)-1;

        weather.G_data =  cellfun(@str2num,table2array(G_table(1+idx_start:1+idx_end, 2:end)));
        weather.fW_data = cellfun(@str2num,table2array(fW_table(1+idx_start:1+idx_end, 2:end)));
        weather.sW_data = cellfun(@str2num,table2array(sW_table(1+idx_start:1+idx_end, 2:end)));
        weather.temp_data = cellfun(@str2num,table2array(temp_table(1+idx_start:1+idx_end, 2:end)));
        weather.rho_data = cellfun(@str2num,table2array(rho_table(1+idx_start:1+idx_end, 2:end)));
        weather.timeline = weather.timeline(idx_start:idx_end);
        weather.timeline_string = weather.timeline_string(idx_start:idx_end,:);
        weather.fulldate_string = weather.fulldate_string(idx_start:idx_end,:);
    end

    
    % Get data after your actual time (remove previous data) from the floor time wrt to yours; e.g. if you are at 9.35, get data from 9.30
    idx_overtime = find(weather.timeline > t_0-7.5*60);
    weather.timeline = weather.timeline(idx_overtime,:);
    weather.timeline_string = weather.timeline_string(idx_overtime:end,:);
    weather.fulldate_string = weather.fulldate_string(idx_overtime:end,:);
    weather.G_data = weather.G_data(idx_overtime,:);
    weather.fW_data = weather.fW_data(idx_overtime,:);
    weather.sW_data = weather.sW_data(idx_overtime,:);
    weather.rho_data = weather.rho_data(idx_overtime,:);
    weather.temp_data = weather.temp_data(idx_overtime,:);

    % Extract the year, month, and day components
    dateStr = weather.fulldate_string(1, 1:10);
    
    % Convert the date string to a date object
    dateObj = datetime(dateStr, 'InputFormat', 'yyyy-MM-dd');
    
    % Format the date object as 'dd-MM-yyyy' and store it in the cell array
    weather.ddmmyy_string = datestr(dateObj, 'dd-mm-yyyy');

    %% Clear temporary variables
    clear opts
end