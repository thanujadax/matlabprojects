function [y,errBar] = getMeanInterpolationCurves(matFilePath,fileStr,calibrationInds)

% Reads the .mat files saved in matFilePath. Each of those files should
% correspond to different interpolation curves
% e.g. curve1 = mean(xcorrMat_1,1);
% errBar is the standard deviation of each curve
% file names should be of pattern fileStr%d.mat

% matFilePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s108';
% fileStr = 'xcorrMat';

% calibrationInds - vector containing the indices of the calibration methods to be used. If empty, all are used.

matFileNames = '*.mat';
matFileNames = fullfile(matFilePath,matFileNames);

matFileDir = dir(matFileNames);

numCurves = length(matFileDir);

if(numCurves==0)
    error('No .mat files found in %s',matFilePath);
end

if(isempty(calibrationInds))
    calibrationInds = 1:numCurves;
end

% load mat file to get the number of data points
load(fullfile(matFilePath,matFileDir(1).name));
numDataPoints = size(xcorrMat,2);

y = zeros(numCurves,numDataPoints);
errBar = zeros(numCurves,numDataPoints);

for i=1:numel(calibrationInds)
    load(fullfile(matFilePath,matFileDir(calibrationInds(i)).name));
    y(i,:) = mean(xcorrMat,1);
    errBar(i,:) = std(xcorrMat);
    clear xcorrMat
end
