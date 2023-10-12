%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 08.09.2023
%% Description 
% This function get the latest weather data in the folder 'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\Forecast\'
% INPUT: 
% "folderDirectory": string, should always be 'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\Forecast\'

% OUTPUT : 
% "latestFolderpath" : string, returns the folder with the latest weather
% forecast data to be used in "load_weather"
%%
function latestFolderpath = get_latest_weather(folderDirectory)
    % List all folders in the directory
    folders = dir(fullfile(folderDirectory, '*_*_*'));
    
    % Extract date and time information from folder names using regular expressions
    folderNames = {folders.name};
    dateTimes = regexp(folderNames, '(\d{6})_(\d{6})', 'tokens');
    
    % Convert date and time strings to numerical values
    dateTimes = cellfun(@(x) str2double([x{1}{1}, x{1}{2}]), dateTimes);
    
    % Find the index of the folder with the latest date and time
    [~, latestIndex] = max(dateTimes);
    
    % Get the name of the latest folder
    latestFolder = folders(latestIndex).name;
    
    latestFolderpath = [folderDirectory,latestFolder,'\preprocess\'];
    % Now, latestFolder contains the name of the folder with the latest date and time
    fprintf('The latest folder is: %s\n', latestFolder);    
end