function [y,errBar] = getMeanInterpolationCurves(matFilePath,fileStr)

% Reads the .mat files saved in matFilePath. Each of those files should
% correspond to different interpolation curves
% e.g. curve1 = mean(xcorrMat_1,1);
% errBar is the standard deviation of each curve
% file names should be of pattern fileStr%d.mat

% matFilePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s108';
% fileStr = 'xcorrMat';

matFileNames = '*.mat';
matFileNames = fullfile(matFilePath,matFileNames);

matFileDir = dir(matFileNames);

numCurves = length(matFileDir);

if(numCurves==0)
    error('No .mat files found in %s',matFilePath);
end

% load mat file to get the number of data points
load(fullfile(matFilePath,matFileDir(1).name));
numDataPoints = size(xcorrMat,2);

y = zeros(numCurves,numDataPoints);
errBar = zeros(numCurves,numDataPoints);

for i=1:numCurves
    load(fullfile(matFilePath,matFileDir(i).name));
    y(i,:) = mean(xcorrMat,1);
    errBar(i,:) = std(xcorrMat);
    clear xcorrMat
end
