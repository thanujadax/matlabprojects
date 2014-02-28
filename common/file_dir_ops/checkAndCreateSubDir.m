function checkAndCreateSubDir(pathToDir,subDir)

% check if subDir exists in pathToDir

% create subDir if it doesn't already exist

% if the last character of subDir is '/', remove it
if(subDir(end)=='/')
    subDir(end) = [];
end

fPath = fullfile(pathToDir,subDir);

if(~isequal(exist(fPath, 'dir'),7)) 
    % 7 = directory
    % doesn't already exist -> create
    status = mkdir(pathToDir,subDir);
    
    if(status)
        % successfully created
        msg = strcat('Subdirectory ', subDir,'/ created inside ',pathToDir);
        disp(msg)
    else
        % folder couldn't be created
        msg = strcat('ERROR: subdirectory ', subDir,...
                '/ cannot be created in ',pathToDir);
        disp(msg)
    end
end