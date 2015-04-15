function thicknessEstimates = doThicknessEstimation(...
    calibrationMethod,inputImageStackFileName,outputSavePath)

% do thickness estimation based on one of the following methods to
% calibrate the similarity curve

calibrationMethod = 3;
% 1 - correlation coefficient across ZY sections, along x axis
% 2 - correlation coefficient across XY sections, along x axis
% 3 - SD of XY per pixel intensity difference
% TODO: methods robust against registration problems

%% Parameters
xyResolution = 5; % nm
maxShift = 20;
maxNumImages = 10; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations

inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/cubes/s108_1-200.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/zyCalibration/s108_1-200_20150409';

%% 

if(calibrationMethod == 1)
    disp('Performing xcorr using ZY stack ...')
    xcorrMat = getXcorrZYstack(inputImageStackFileName,maxShift,maxNumImages);
    disp('done!')
    % each column of xcorrMat corresponds to a sequence of shifted frames of
    % the zy plane along the x axis.
    % each row corresponds to one starting image (zy section) of the stack

elseif(calibrationMethod == 2)
    disp('Performing xcorr using shifted XY images stack ...')
    xcorrMat = getXcorrXYstack(inputImageStackFileName,maxShift,maxNumImages);
    disp('done!')    
    
elseif(calibrationMethod == 3)
    disp('Calculating SD of intensity deviation using shifted XY sections stack ...')
    xcorrMat = getIntensityDeviationXYstack(inputImageStackFileName,maxShift,maxNumImages);
    disp('done!')
else
    error('Unrecognized calibration method specifed. Check calibrationMethod')
end

% save output
disp('Saving xcorrMat')
xcorrFileName = fullfile(outputSavePath,'xcorrMat.mat');
save(xcorrFileName,'xcorrMat')
disp('done')

% %% plot
% shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g');

%% Predict
% predict section thickness for the data set
relZresolution = predictThicknessFromCurve(...
        inputImageStackFileName,xcorrMat,maxShift,calibrationMethod);
thicknessEstimates = relZresolution .* xyResolution;
figure;plot(thicknessEstimates)
title('Section thickness estimates (nm)')
xlabel('Section index')
ylabel('Estimated thickness (nm)')
