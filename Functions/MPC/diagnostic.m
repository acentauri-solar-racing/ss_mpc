%% Ngo Tony
% This code was written with MATLAB R2022b. Errors may occur with other
% versions, last updated: 06.09.2023
%% Description 
% This function output wheter there have been any infeasibilities in the
% solver

% INPUT: 
% "max_iter": maximum number of iteration of the mpc simulation
% "dist": distance vector 

% OUTPUT : 
% "diagnostics": diagnostics results

%%

function [diagnostics] = diagnostic(max_iter,dist)
    %DIAGNOSTIC Summary of this function goes here
    %   Detailed explanation goes here
    feasibility         = zeros(max_iter, 1);
    error               = struct;
    flag = 0;
    
    for k = 1:max_iter
        text            = fileread(['diagnostic/diary_' num2str(k)]);
        feas_tmp        = strfind(text, 'EXIT: Optimal Solution Found.');
        if isempty(feas_tmp)
            flag = 1;
            feasibility(k)     = 0;
            tmp                = strfind(text, 'EXIT');
            error.message{k}   = text(tmp:tmp+40);
            error.timestep(k)  = k;
        else
            feasibility(k)     = 1;
        end
        if k == max_iter*dist
            feasibility(k)     = 1;
        end
        % Delete file
        name_tmp        = ['diagnostic/diary_' num2str(k)];
        delete(name_tmp)
    end
    
    notFeasible                 = find(~feasibility);
    diagnostics.alwaysFeasible    = all(feasibility);
    diagnostics.noFeasibility(1,:) = notFeasible;
    diagnostics.noFeasibility(2,:) = dist(notFeasible);
    
    if flag == 1
        diagnostics.errors.message    = error.message(notFeasible); 
        diagnostics.errors.timestep   = error.timestep(notFeasible);
    end

    myFolder = 'diagnostic';
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(myFolder, '*.fig'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(myFolder, baseFileName);
        fprintf(1, 'Now deleting %s\n', fullFileName);
        delete(fullFileName);
    end
end

