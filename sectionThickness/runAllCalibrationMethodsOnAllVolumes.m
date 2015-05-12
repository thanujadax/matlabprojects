function runAllCalibrationMethodsOnAllVolumes...
    (inputImageStackFileName,outputSavePath,params)

% run all calibration methods for one volume and save the calibration
% curves and the predictions in the outputPath

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

imageCubeDirectory = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes';
outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/20150512';

% inputImageStackFileName = '/home/thanuja/projects/data/FIBSEM_dataset/largercubes/s108/s108.tif';
% outputSavePath = '/home/thanuja/projects/tests/thickness/similarityCurves/20150512';

% get the list of directories
[sampleDirectories,~] = subdir(imageCubeDirectory);

% get the tiff stack from each directory. some may have multiple tiff
% stacks.

for i=1:length(sampleDirectories)
    
    sampleSubDirName = sampleDirectories{i};
    % read all image stacks in this sample
    imageStackFileString = strcat('*.',params.imgStackFileExt);
    imageStackFileString = fullfile(sampleSubDirName,imageStackFileString);
    imageStackDir = dir(imageStackFileString);
    
%     str1 = sprintf('Processing image stack %s',sampleSubDirName,calibrationMethod);
%     disp(str1)    
    
    % process each image stack in the sample
    for j=1:length(imageStackDir)
        inputImageStackFileName = fullfile...
            (sampleSubDirName,imageStackDir(j).name);
        
        tokenizedFName = strsplit(inputImageStackFileName,filesep);
        nameOfStack = strtok(tokenizedFName(end),'.');
        nameOfStack = nameOfStack{1};
        % check if subdir exists. if not create.
        checkAndCreateSubDir(outputSavePath,nameOfStack);
        outputSavePath_i = fullfile(outputSavePath,nameOfStack);

        % writes output to output path as txt file. Col vector.
        for calibrationMethod=1:7
            str1 = sprintf('Running calibration method %02d on image stack %s',calibrationMethod,sampleSubDirName);
            disp(str1)
            thicknessEstimates = doThicknessEstimation(...
            calibrationMethod,inputImageStackFileName,outputSavePath_i,params);
        end    
    end    
end

% save parameters for future reference
outputParamsFileName = fullfile(outputSavePath,'params.mat');
save(outputParamsFileName,params);
