function fmOut = readFmFile(pathToSliceFeatures,fileIndex)

% Inputs:
%   fileIndex - sequence number of the file required

% reads the feature matrix for the ith section saved in the directory given
% by pathToSliceFeatures

inputSections = dir(pathToSliceFeatures);

inputFm_FilePath = fullfile(pathToSliceFeatures,inputSections(fileIndex).name);

load(inputFm_FilePath);

fmOut = fm;