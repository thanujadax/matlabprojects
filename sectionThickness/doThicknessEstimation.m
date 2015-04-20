function thicknessEstimates = doThicknessEstimation(...
    calibrationMethod,inputImageStackFileName,outputSavePath)

% do thickness estimation based on one of the following methods to
% calibrate the similarity curve

calibrationMethod = 5;
% 1 - correlation coefficient across ZY sections, along x axis
% 2 - correlation coefficient across XY sections, along x axis
% 3 - SD of XY per pixel intensity difference
% TODO: methods robust against registration problems

%% Parameters
xyResolution = 5; % nm
maxShift = 20;
maxNumImages = 20; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
numPairs = 2; % number of section pairs to be used to estimate the thickness of one section

inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/cubes/s108_1-200.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s108_1-200_20150416';

%% 

if(calibrationMethod == 1)
    disp('Calculating c.o.c decay curve using ZY stack, along X ...')
    xcorrMat = getXcorrZYstack(inputImageStackFileName,maxShift,maxNumImages);
    disp('done!')
    % each column of xcorrMat corresponds to a sequence of shifted frames of
    % the zy plane along the x axis.
    % each row corresponds to one starting image (zy section) of the stack

elseif(calibrationMethod == 2)
    disp('Calculating c.o.c decay curve using XY images stack, along Y ...')
    xcorrMat = getXcorrXYstack(inputImageStackFileName,maxShift,maxNumImages);
    disp('done!')    
    
elseif(calibrationMethod == 3)
    disp('Calculating SD of intensity deviation curve using shifted XY sections stack, along X ...')
    xcorrMat = getIntensityDeviationXYstack(inputImageStackFileName,maxShift,maxNumImages);
    disp('done!')
elseif(calibrationMethod == 4)
    disp('Calculating c.o.c decay curve using ZY stack, along Y ...')
    xcorrMat = getXcorrZYstackY(inputImageStackFileName,maxShift,maxNumImages);
    disp('curve estimation done')
elseif(calibrationMethod == 5)
    disp('Calculating c.o.c decay curve using XY images stack, along X ...')
    xcorrMat = getXcorrXYstackX(inputImageStackFileName,maxShift,maxNumImages);
    disp('curve estimation done!')    
    
elseif(calibrationMethod == 6)
    disp('Calculating c.o.c decay curve using XZ images stack, along X ...')
    xcorrMat = getXcorrXZstackX(inputImageStackFileName,maxShift,maxNumImages);
    disp('curve estimation done!')    
elseif(calibrationMethod == 7)
    disp('Calculating c.o.c decay curve using XZ images stack, along Y ...')
    xcorrMat = getXcorrXZstackY(inputImageStackFileName,maxShift,maxNumImages);
    disp('curve estimation done!')
else
    error('Unrecognized calibration method specified. Check calibrationMethod')
end

% save output
disp('Saving xcorrMat')
matName = sprintf('xcorrMat%d.mat',calibrationMethod);
xcorrFileName = fullfile(outputSavePath,matName);
save(xcorrFileName,'xcorrMat')
disp('done')

% %% plot
% shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g');

%% Predict
% predict section thickness for the data set
% relZresolution = predictThicknessFromCurve(...
%         inputImageStackFileName,xcorrMat,maxShift,calibrationMethod);
    
relZresolution = predictThicknessFromCurveFromMultiplePairs(...
        inputImageStackFileName,xcorrMat,maxShift,calibrationMethod,numPairs);
% each row contains one set of estimates for all the sections. First row is
% from the images i and i+1, the second row is from images i and i+2 etc    
thicknessEstimates = relZresolution .* xyResolution;
if(size(thicknessEstimates,1)==1)
    figure;plot(thicknessEstimates);
elseif(size(thicknessEstimates,1)==2)
    figure;plot(thicknessEstimates(1,:),'r');
    hold on
    plot((thicknessEstimates(2,:).* 0.5),'g');
    hold off
end
title('Section thickness estimates (nm)')
xlabel('Section index')
ylabel('Estimated thickness (nm)')
