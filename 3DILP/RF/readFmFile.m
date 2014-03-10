function fmOut = readFmFile(pathToSliceFeatures,fileIndex,fileName)

% Inputs:
%   fileIndex - sequence number of the file required
%   fileName - if the specific file name is know
% reads the feature matrix for the ith section saved in the directory given
% by pathToSliceFeatures. fileIndex is irrelevant if fileName is set.

% if file name is not given, it loads the file given by fileIndex,
% considering all the files available in the directory.

fileNameStr = '/*.mat';

if(nargin<3)
    fileName = [];
end

if(isempty(fileName))

    inputSections = dir(strcat(pathToSliceFeatures,fileNameStr));

    inputFm_FilePath = fullfile(pathToSliceFeatures,inputSections(fileIndex).name);

    fmOut = importdata(inputFm_FilePath);


    
else
    % read the file given by file name
    fileNameToRead = fullfile(pathToSliceFeatures,fileName);
    fmOut = importdata(fileNameToRead);
end