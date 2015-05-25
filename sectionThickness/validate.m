function validate()

%% Inputs
inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes/s108/s108.tif';
precomputedMatFilePath = '/home/thanuja/projects/tests/thickness/similarityCurves/20150525/s108';
outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/20150525/s108/validation';
fileStr = 'xcorrMat'; % general string that defines the .mat file
distMin = 0;
saveOnly = 0;
xResolution = 5; % nm
yResolution = 5; % nm

validateUsingXresolution = 1 ; % if 0, validation is done using y resolution.
% For FIBSEM, the resolution of x and y are known (~ 5 nm)
% if set we treat y as the known resolution to calibrate the decay curves
% and x as the direction for which  

%% Param
method = 'spline'; % method of interpolation
method2 = 'linear';

color = 'g';
predictionFigureFileStr = 'prediction';

%% Validation

if(validateUsingXresolution)
    % predict y resolution using x
    calibrationMethods = [1 3 5];
    calibrationString = 'Avg c.o.c decay using X resolution';
    calibrationFigureFileString = 'coc_xResolution_ensemble';
    subTitle = 'XZ_y';
else
    % predict x resolution using y
    calibrationMethods = [2 4 6];
    calibrationString = 'Avg c.o.c decay using Y resolution';
    calibrationFigureFileString = 'coc_yResolution_ensemble';
    subTitle = 'YZ_x';
end

% get the calibration curves from the precomputed .mat files        
% get the avg calibration curve
[meanVector,stdVector] = makeEnsembleDecayCurveForVolume...
    (precomputedMatFilePath,fileStr,0,calibrationMethods);

% plot decay curve
plotSaveMeanCalibrationCurveWithSD...
    (inputImageStackFileName,calibrationString,saveOnly,...
    distMin,meanVector,stdVector,color,outputSavePath,calibrationFigureFileString);

%% method1 predict the thickness/resolution in the assumed unknown direction
if(validateUsingXresolution)
    predictedThickness = predictThicknessXZ_Y...
            (inputImageStackFileName,meanVector,xResolution,...
            distMin,method);
    
else
    predictedThickness = predictThicknessYZ_X...
            (inputImageStackFileName,meanVector,xResolution,...
            distMin,method);
end
% plot predicted thickness
plot(predictedThickness);
titleStr = sprintf('Predicted thickness %s (%s interploation)',subTitle,method);
title(titleStr)
xlabel('Inter-section interval')
ylabel('Thickness (nm)')
% save plot
predictionFileName = sprintf('%s_%s_%s',predictionFigureFileStr,subTitle,method);
predictionFileName = fullfile(outputSavePath,predictionFileName);
print(predictionFileName,'-dpng');
% save thickness in txt file
save(fullfile(predictionFileName,'.dat'),'predictedThickness','-ASCII');
% calculate the error, mean error and the variance

%% method 2: predict the thickness/resolution in the assumed unknown direction
predictedThickness = predictThicknessYZ_X...
        (inputImageStackFileName,meanVector,xResolution,...
        distMin,method2);

% plot predicted thickness
plot(predictedThickness);
titleStr = sprintf('Predicted thickness %s (%s interploation)',subTitle,method2);
title(titleStr)
xlabel('Inter-section interval')
ylabel('Thickness (nm)')
% save plot
predictionFileName = sprintf('%s_%s_%s',predictionFigureFileStr,subTitle,method2);
predictionFileName = fullfile(outputSavePath,predictionFileName);
print(predictionFileName,'-dpng');
% save thickness in txt file
save(fullfile(predictionFileName,'.dat'),'predictedThickness','-ASCII');
% calculate the error, mean error and the variance

