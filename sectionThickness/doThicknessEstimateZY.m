function thicknessEstimates = doThicknessEstimateZY()
% produce similarity decay curve using zy sections to predict thickness of
% xy sections

%% Parameters
xyResolution = 5; % nm
maxShift = 20;
maxNumImages = 10;

% inputImageStackFileName = '/home/thanuja/projects/data/MSB/substacks/A_smaller_1-250.tif';
% outputSavePath = '/home/thanuja/projects/tests/thickness/zyCalibration/test20150407';
inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/cubes/s108_1-200.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/zyCalibration/s108_1-200_20150409';

%% 
disp('Performing xcorr using ZY stack ...')
xcorrMat = getXcorrZYstack(inputImageStackFileName,maxShift,maxNumImages);
disp('done!')
% each column of xcorrMat corresponds to a sequence of shifted frames of
% the zy plane along the x axis.
% each row corresponds to one starting image (zy section) of the stack

% save output
disp('Saving xcorrMat')
xcorrFileName = fullfile(outputSavePath,'xcorrMat.mat');
save(xcorrFileName,'xcorrMat')
disp('done')

%% plot
shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g');

%% Predict
% predict section thickness for the data set
relZresolution = predictThicknessFromCurve(...
        inputImageStackFileName,xcorrMat,maxShift);
thicknessEstimates = relZresolution .* xyResolution;

