% thickness estimation batch srcript
function main_estimateThicknessBatch()

%% Parameters
calibrationMethod = 5;
% 1 - correlation coefficient across ZY sections, along x axis
% 2 - correlation coefficient across XY sections, along x axis
% 3 - SD of XY per pixel intensity difference
% 4 - c.o.c across ZY sections along Y
% 5 - c.o.c across XY sections, along X
% 6 - c.o.c acroxx XZ sections, along X
% 7 - 
% TODO: methods robust against registration problems
param.imgStackFileExt = 'tif';
params.xyResolution = 5; % nm
params.maxShift = 15;
params.maxNumImages = 50; % number of sections to initiate calibration.
                % the calibration curve is the mean value obtained by all
                % these initiations
params.numPairs = 2; % number of section pairs to be used to estimate the thickness of one section
params.plotOutput = 0; % don't plot intermediate curves.
params.usePrecomputedCurve = 0;
params.pathToPrecomputedCurve = '';

% inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/cubes/s108_1-200.tif';
outputSavePath = '/home/thanuja/projects/tests/thickness/batchEstimation/20150504';

%% File paths
imageCubeDirectory = '/home/thanuja/projects/data/FIBSEM_dataset/cubes';

% get the list of directories
[sampleDirectories,~] = subdir(imageCubeDirectory);

% tiff stacks are storedd in each directory
for i=1:length(sampleDirectories)
    sampleSubDirName = sampleDirectories{i};
    % read all image stacks in this sample
    imageStackFileString = strcat('*.',param.imgStackFileExt);
    imageStackFileString = fullfile(sampleSubDirName,imageStackFileString);
    imageStackDir = dir(imageStackFileString);
    
    % process each image stack in the sample
    for i=1:length(imageStackDir)
        inputImageStackFileName = fullfile...
            (sampleSubDirName,imageStackDir(i).name);
        thicknessEstimates = doThicknessEstimation(...
    calibrationMethod,inputImageStackFileName,outputSavePath,params);
    
    end
    
    % TODO: do we use the same curve for ths same sample?
    % makes sense to do so. Make aggregate curve?
end