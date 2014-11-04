% Vesicle filter 3D

%% Parameters
innerRadius = 1;
outerRadius = 3;

%% File paths
% tif file
inputImageVolumeFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/GAS1_mss1_small_50.tif';
outputImageVolumeFileName = '/home/thanuja/Dropbox/PROJECTS/vesicleFilter/output/small50_1.tif';

%% Read file
disp('Reading input image file')
disp('inputImageVolumeFileName');
inputImageStack = readTiffStackToArray(inputImageVolumeFileName);
inputImageStack = inputImageStack./255;
% % test file read
% writeImageMatrixToTiffStack(inputImageStack,'inputStack.tiff');
%% convolution with 3D structure element (vesicle)

disp('Generating vesicle structure element')
vesicleElement3D = getVesicleElement(innerRadius,outerRadius);

% % test vesicle element
% outputFileName = 'vesicleElement.tiff';
% writeImageMatrixToTiffStack(vesicleElement3D,outputFileName);

% convn (or fftn for fft based convolution)
disp('convolution...')
convolutionResult = convn(inputImageStack,vesicleElement3D,'same');
% interpn: to interpolate the convolution result at the points you want

%% visualize
disp('Writing output to')
disp(outputImageVolumeFileName);
writeImageMatrixToTiffStack(convolutionResult,outputImageVolumeFileName);
disp('Done :-)')
% currently to visualize the tiff volume, we are using fiji