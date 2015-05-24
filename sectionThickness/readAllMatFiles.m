function combinedMat = readAllMatFiles(matFilePath,fileStr,z)

% read each matfile and combine into one
% what's the mat file structure?
%   each row corresponds to a different image
%   each column corresponds to a different distance
%   so all the mat files can be appended beneath the previous one

matFileNames = strcat(fileStr,'*.mat');
matFileNames = fullfile(matFilePath,matFileNames);

matFileDir = dir(matFileNames);
TotNumCurves = length(matFileDir);

if(z)
    % along z axis. only 3 curves to be combined
    indices = [8,9,10];
    % numCurves = numel(indices); % 8,9,10
    disp('Reading .mat files ending with 08, 09 and 10 - along z axis')
    
else
    % x,y directions. i.e. 6 curves
    indices = [1,2,4,5,6,7];
    % numCurves = 6; % 1,2,4,5,6,7
    disp('Reading .mat files ending with 1,2,4,5,6,7 - along x,y')
        
end

combinedMat = readMats(indices,matFilePath,matFileDir);


function mat = readMats(indices,matFilePath,matFileDir)

mat = [];
numInds = length(indices);

for i=1:numInds
    clear xcorrMat
    load(fullfile(matFilePath,matFileDir(indices(i)).name));
    mat = [mat; xcorrMat];
    
end