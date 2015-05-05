function runAllCalibrationMethodsOnVolume...
    (inputImageStackFileName,outputSavePath,params)

% run all calibration methods for one volume and save the calibration
% curves and the predictions in the outputPath

params.predict = 1; % set to 0 if only the interpolation curve is required.
params.xyResolution = 5; % nm
params.maxShift = 15;
params.maxNumImages = 60; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
params.numPairs = 1; % number of section pairs to be used to estimate the thickness of onesection
params.plotOutput = 0;
params.usePrecomputedCurve = 0;
params.pathToPrecomputedCurve = '';

inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes/s108/s108.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s108';

for calibrationMethod=1:7
    str1 = sprintf('Running calibration method %02d',calibrationMethod);
    disp(str1)
    thicknessEstimates = doThicknessEstimation(...
    calibrationMethod,inputImageStackFileName,outputSavePath,params);
    
end