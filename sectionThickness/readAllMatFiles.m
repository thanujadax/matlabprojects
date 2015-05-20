function combinedMat = readAllMatFiles(matFilePath,fileStr,z)

% read each matfile and combine into one
% what's the mat file structure?
%   each row corresponds to a different image
%   each column corresponds to a different distance
%   so all the mat files can be appended beneath the previous one

matFileNames = '*.mat';
matFileNames = fullfile(matFilePath,matFileNames);

matFileDir = dir(matFileNames);
TotNumCurves = length(matFileDir);

if(z)
    % along z axis. only 3 curves to be combined
    numCurves = 3; % 8,9,10
else
    % x,y directions. i.e. 6 curves
    numCurves = 6; % 1,2,4,5,6,7
end
