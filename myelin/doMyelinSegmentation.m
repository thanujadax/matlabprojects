function doMyelinSegmentation(inputFilePaths,params,outputPath)

% Inputs
%   raw images
%   label images (for error calculation)

% Outputs
%   save myelin segmentation
%       1. as separate label images
%       2. as overlay on raw images
%   error measures

%% File paths
inputFilePaths.fileExt = 'png';
inputFilePaths.rawImages = '/home/thanuja/projects/data/myelin/ssSEM/s909/raw20150428';
inputFilePaths.myelinProbabilityMaps = '/home/thanuja/projects/data/myelin/ssSEM/s909/probabilityMaps20150428';
outputPath = '/home/thanuja/projects/data/myelin/ssSEM/s909/output20150428';

%% Parameters
params.sigmaGaussianBlur = 2;
params.maskSizeGaussianBlur = 7;
params.threshold_pixelIntensity = 0.15; % to generate binary images for myelin labels
params.numPixelsCC = 3000; % minimum number of pixels for CC to be considered as myelin

%% Perform segmentation on raw images
inputMylProbMapFilesString = strcat('*.',inputFilePaths.fileExt);
inputMylProbMapFilesString = fullfile(inputFilePaths.myelinProbabilityMaps,...
                                    inputMylProbMapFilesString);
inputProbabilityMapDirectory = dir(inputMylProbMapFilesString);

for i = 1:length(inputProbabilityMapDirectory)
    probabilityMapFileName = fullfile(inputFilePaths.myelinProbabilityMaps,inputProbabilityMapDirectory(i).name);
    str1 = strcat('Reading probability map ',probabilityMapFileName);
    disp(str1)
    probabilityMap_i = imread(probabilityMapFileName);
    %figure; imshow(probabilityMap_i);
    disp('Segmenting myelin...')
    myelinSegments = getMyelinSegments1(probabilityMap_i,params);
    disp('done!')
    % save myelin segmentation
    outputFileName = sprintf('segmentedMyelin%d.png',i);
    outputFileName = fullfile(outputPath,outputFileName);
    str1=strcat('Saving output to file: ',outputFileName);
    disp(str1);
    imwrite(myelinSegments,outputFileName);
    disp('done');
end