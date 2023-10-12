%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 08.09.2023
%% Description 
% This function returns the polinomial fit parameters for the solar
% irradiation G, assuming it's only depending on the time

% INPUT: 

% OUTPUT : 
% t_0: time in seconds of the actual time

%%
function t_0 = get_machine_time_s()
    format = 'HHMM';
    hourstr = datestr(datetime,format);
    t_0 = str2num(hourstr(1:2))*60*60+str2num(hourstr(3:4))*60;
end