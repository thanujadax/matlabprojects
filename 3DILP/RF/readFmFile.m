function fmOut = readFmFile(pathToSliceFeatures,i)

% reads the feature matrix for the ith section saved in the directory given
% by pathToSliceFeatures

inputSections = dir(pathToSliceFeatures);

inputFm_FilePath = fullfile(pathToSliceFeatures,inputSections(i).name);

load(inputFm_FilePath);

fmOut = fm;