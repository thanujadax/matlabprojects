% thickness estimation batch srcript
function main_estimateThicknessBatch()

%% Parameters
calibrationMethod = 1;
% % 1 - c.o.c across ZY sections, along x axis
% % 2 - c.o.c across XY sections, along Y axis
% % 3 - SD of XY per pixel intensity difference
% % 4 - c.o.c across ZY sections along Y
% % 5 - c.o.c across XY sections, along X
% % 6 - c.o.c acroxx XZ sections, along X
% % 7 - c.o.c acroxx XZ sections, along Y
% TODO: methods robust against registration problems
param.imgStackFileExt = 'tif';

params.predict = 1;
params.xyResolution = 5; % nm
params.maxShift = 10;
params.maxNumImages = 100; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
params.numPairs = 1; % number of section pairs to be used to estimate the thickness of one section
params.plotOutput = 0; % don't plot intermediate curves.
params.usePrecomputedCurve = 0;
params.pathToPrecomputedCurve = '';

%% File paths

outputSavePath = '/home/thanuja/projects/tests/thickness/batchEstimation/20150504/calibration01';
imageCubeDirectory = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes';

% get the list of directories
[sampleDirectories,~] = subdir(imageCubeDirectory);

% tiff stacks are stored in each directory
parfor i=1:length(sampleDirectories)
    
    sampleSubDirName = sampleDirectories{i};
    % read all image stacks in this sample
    imageStackFileString = strcat('*.',param.imgStackFileExt);
    imageStackFileString = fullfile(sampleSubDirName,imageStackFileString);
    imageStackDir = dir(imageStackFileString);
    
    str1 = sprintf('Processing image stack %s with method %d',sampleSubDirName,calibrationMethod);
    disp(str1)
    % process each image stack in the sample
    for j=1:length(imageStackDir)
        inputImageStackFileName = fullfile...
            (sampleSubDirName,imageStackDir(j).name);
        thicknessEstimates = doThicknessEstimation(...
    calibrationMethod,inputImageStackFileName,outputSavePath,params);
    % writes output to output path as txt file. Col vector.
    end
    
    % TODO: do we use the same curve for ths same sample?
    % makes sense to do so. Make aggregate curve?
end