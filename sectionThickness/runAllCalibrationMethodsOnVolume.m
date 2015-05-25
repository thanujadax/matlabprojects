function runAllCalibrationMethodsOnVolume...
    (inputImageStackFileName,outputSavePath,params)

% run all calibration methods for one volume and save the calibration
% curves and the predictions in the outputPath

% % calibrationMethod
% % 1 - c.o.c across XY sections, along X
% % 2 - c.o.c across XY sections, along Y axis
% % 3 - c.o.c across ZY sections, along x axis
% % 4 - c.o.c across ZY sections along Y
% % 5 - c.o.c acroxx XZ sections, along X
% % 6 - c.o.c acroxx XZ sections, along Y
% % 7 - c.o.c across XY sections, along Z
% % 8 - c.o.c across ZY sections, along Z
% % 9 - c.o.c. across XZ sections, along Z
% % 10 - SD of XY per pixel intensity difference
params.predict = 1; % set to 0 if only the interpolation curve is required.
params.xyResolution = 5; % nm
params.maxShift = 10;
params.minShift = 0;
params.maxNumImages = 700; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
params.numPairs = 1; % number of section pairs to be used to estimate the thickness of onesection
params.plotOutput = 1;
params.suppressPlots = 1;
params.usePrecomputedCurve = 0;
params.pathToPrecomputedCurve = '';
params.imgStackFileExt = 'tif';

inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes/s108/s108.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/20150525/s108';

for calibrationMethod=1:10
    str1 = sprintf('Running calibration method %02d',calibrationMethod);
    disp(str1)
    thicknessEstimates = doThicknessEstimation(...
    calibrationMethod,inputImageStackFileName,outputSavePath,params);
    
end