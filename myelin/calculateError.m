function errors = calculateError()

% Calculates errors of the segmentations
%   errors.GCE - Globabl Consistency Error (Martin et al 2001)
%   errors.LCE - Local Consistency Error (Martin et al 2001)
%   errors.F1score - 

%% File names and paths
inputFilePaths.groundTruthDirectory = '/home/thanuja/projects/data/myelin/ssSEM/s909/groundTruth20150428';
inputFilePaths.segmentationDirectory = '/home/thanuja/projects/data/myelin/ssSEM/s909/output20150428';
inputFilePaths.inputFileType = 'png';
outputFileName = '/home/thanuja/projects/data/myelin/ssSEM/s909/output20150428/predictionError.txt';

%% Parameters

%% Error calculations
% compare if there is an equal number of files. If not equal, pick the
% least and report as warning.

inputGroundTruthFilesString = strcat('*.',inputFilePaths.inputFileType);
inputGroundTruthFilesString = fullfile...
    ('/home/thanuja/projects/data/myelin/ssSEM/s909/groundTruth20150428',inputGroundTruthFilesString);
inputGroundTruthDirectory = dir(inputGroundTruthFilesString);

inputSegmentationFilesString = strcat('*.',inputFilePaths.inputFileType);
inputSegmentationFilesString = fullfile...
    ('/home/thanuja/projects/data/myelin/ssSEM/s909/output20150428',inputSegmentationFilesString);
inputSegmentationDirectory = dir(inputSegmentationFilesString);

numGroundTruthImages = length(inputGroundTruthDirectory);
numSegmentImages = length(inputSegmentationDirectory);

str1 = sprintf('Num ground truth images found: %d', numGroundTruthImages);
disp(str1);
str1 = sprintf('Num segmented images found: %d', numSegmentImages);
disp(str1);

if(numGroundTruthImages==numSegmentImages)
    numImagesToProcess = numGroundTruthImages;
else
    disp('WARNING: calculateError.m: the number of ground truth images and the number of segmented images do not match.!');
    numImagesToProcess = min(numGroundTruthImages,numSegmentImages);
    str1('Comparing the first %d images.',numImagesToProcess);
    disp(str1)
end
errors.ri = zeros(1,numImagesToProcess);
errors.gce = zeros(1,numImagesToProcess);
errors.vi = zeros(1,numImagesToProcess);

for i=1:numImagesToProcess
    str1 = sprintf('Processing image pair %d out of %d',i,numImagesToProcess);
    disp(str1)
    
    % get ground truth image
    groundTruthImageFileName = fullfile...
        (inputFilePaths.groundTruthDirectory,inputGroundTruthDirectory(i).name);
    groundTruthImage = imread(groundTruthImageFileName);
    if(size(groundTruthImage,3)>1)
        % convert RGB to indexed image
        [groundTruthImage,map] = rgb2ind(groundTruthImage,65536);
    end
    % get segmented image
    segmentedImageFileName = fullfile...
        (inputFilePaths.segmentationDirectory,inputSegmentationDirectory(i).name);
    segmentedImage = imread(segmentedImageFileName);
    if(size(segmentedImage,3)>1)
        % convert RGB to indexed image
        [segmentedImage,map] = rgb2ind(segmentedImage,65536);
    end
    
    % calculate various errors
    [errors.ri(i),errors.gce(i),errors.vi(i)] = compareSegmentations...
                        (segmentedImage,groundTruthImage);
                    
end

%% write errors to file
fileID = fopen(outputFileName,'w');
fprintf(fileID,'ImageID \t RandIndex \t GCE \t\t VoI \n');
for i=1:numImagesToProcess
    fprintf(fileID,'%d \t \t %0.3e \t %0.3e \t %0.3e \n',i,errors.ri(i),errors.gce(i),errors.vi(i));
end
fclose(fileID);
