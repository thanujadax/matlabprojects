function [y,errBar] = makeEnsembleDecayCurveForVolume(matFilePath,fileStr,zDirection)

% Inputs:
% matFilePath = '/home/thanuja/projects/tests/thickness/similarityCurves/20150512/s108';
% fileStr = 'xcorrMat'; % general string that defines the .mat file
% zDirection = 1; % 0 for x,y direction decay curves

if(zDirection)
    directionStr = 'z direction';
    plotColor = 'b';
else
    directionStr = 'x,y directions';
    plotColor = 'g';
end

% read relevant mat files

% combine them
% get mean and variance
% plot smoot curve - using interpolation ?

combinedMat = readAllMatFiles(matFilePath,fileStr,zDirection);
[numR,numC] = size(combinedMat);
% each row corresponds to a different image. each column correspond to c.c.
% at different distances

y = mean(combinedMat,1);
errBar = std(combinedMat);

%% Plot
tokenizedSubDirName = strsplit(matFilePath,filesep);
volName = tokenizedSubDirName{end};
titleStr = sprintf('Ensemble similarity decay curve: vol %s - %s',volName,directionStr);
xlabelStr = 'Distance (#px or #sections)';
ylabelStr = 'Coefficient of Correlation';
transparent = 0;
figure;
shadedErrorBar((0:(numC-1)),y,errBar,plotColor,transparent,...
    titleStr,xlabelStr,ylabelStr);