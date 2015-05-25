function validate()

inputImageStackFileName = '';
precomputedMatFilePath = '';
outputSavePath = '';
fileStr = 'xcorrMat'; % general string that defines the .mat file

validateUsingXresolution = 0 ; % if 0, validation is done using y resolution.
% For FIBSEM, the resolution of x and y are known (~ 5 nm)
% if set we treat y as the known resolution to calibrate the decay curves
% and x as the direction for which  

if(validateUsingXresolution)
    % predict y resolution using x
    calibrationMethods = [1 3 5];
else
    % predict x resolution using y
    calibrationMethods = [2 4 6];
end

% get the calibration curves from the precomputed .mat files        
% get the avg calibration curve
[meanVector,errBar] = makeEnsembleDecayCurveForVolume...
    (precomputedMatFilePath,fileStr,0,calibrationMethods);    

% predict the thickness/resolution in the assumed unknown direction


% calculate the error, mean error and the variance