function thicknessEstimates = doThicknessEstimateZY()
% produce similarity decay curve using zy sections to predict thickness of
% xy sections

%% Parameters
xyResolution = 5; % nm
maxShift = 20;
maxNumImages = 10;
inputImageStackFileName = 'name.tif';
outputSavePath = 'outputPath';

%% 
xcorrMat = getXcorrZYstack(inputImageStackFileName,maxShift,maxNumImages);
% each column of xcorrMat corresponds to a sequence of shifted frames of
% the zy plane along the x axis.
% each row corresponds to one starting image (zy section) of the stack

% save output

%% plot
shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g');

%% Predict
% predict section thickness for the data set
