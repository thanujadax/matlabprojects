function plotSaveMeanCalibrationCurveWithSD...
    (inputImageStackFileName,calibrationString,saveOnly,...
    distMin,meanVector,stdVector,color,outputSavePath,calibrationFigureFileString)

tokenizedFName = strsplit(inputImageStackFileName,filesep);
nameOfStack = strtok(tokenizedFName(end),'.');
nameOfStack = nameOfStack{1};


titleStr = sprintf('Similarity curve: %s. Vol %s',calibrationString,nameOfStack);
title(titleStr)
xlabelStr = 'Shifted pixels';
ylabelStr = 'Coefficient of Correlation';
transparent = 0;
if(saveOnly)
    set(gcf,'Visible','off');
end
distMax = numel(meanVector);
shadedErrorBar((distMin:(distMax-1)),meanVector,stdVector,color,transparent,...
    titleStr,xlabelStr,ylabelStr);
% save calibration curve figure
% calibrationFigureFileName = strcat(calibrationFigureFileString,'.png')
calibrationFigureFileName = fullfile(outputSavePath,calibrationFigureFileString);
print(calibrationFigureFileName,'-dpng');
