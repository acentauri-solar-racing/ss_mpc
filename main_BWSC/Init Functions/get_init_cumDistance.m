function cumDistance = get_init_cumDistance()
    % 1. Get the current date in the format 'YYYYMMDD'
    currentDate = datestr(now, 'yyyymmdd');
    
    % 2. Create the file name with the prefix 'GPS'
    fileName = strcat(currentDate, '_GPS');
    
    % 3. Specify the directory path
    directoryPath = 'G:\Shared drives\AlphaCentauri\SolarCar_22 23\6. Strategy & Simulation\ss_online_data\Solar_car\GPS';
    
    % 4. Combine the directory path and the file name
    fullFilePath = fullfile(directoryPath, fileName);
    
    % 5. Read the CSV file
    data = readtable(fullFilePath, 'Delimiter', ',');

    % 6. Getting latest value
    cumDistance = data(end,4);
end
