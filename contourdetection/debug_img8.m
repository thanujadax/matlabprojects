% debug image 08

rawType = 'tif';
neuronProbabilityType = 'png';
membraneProbabilityType = 'tiff';
mitoProbabilityType = 'png';

inputPath = '/home/thanuja/projects/data/toyData/set8/';
outputPath = '/home/thanuja/projects/tests/contours/20150514';
outputPathPNG = '/home/thanuja/projects/tests/contours/20150514png';
% read all images in the raw images file path
rawImagePath = fullfile(inputPath,'raw');
allRawFiles = dir(fullfile(rawImagePath,'*.tif'));

% for each file
numFiles = length(allRawFiles);
%for i=1:numFiles
i = 9;
    disp(i);
    imageFileName = allRawFiles(i).name;
    segmentationOut = doILP_w_dir(inputPath,imageFileName,i,...
        rawType,neuronProbabilityType,membraneProbabilityType,mitoProbabilityType);
    % save segmentation output
    writeFileName = fullfile(outputPath,imageFileName);
    imwrite(segmentationOut,writeFileName,'tif');
    pngFileName = sprintf('%d.png',(i-1));
    pngFileName = fullfile(outputPathPNG,pngFileName);
    imwrite(segmentationOut,pngFileName,'png');
%end