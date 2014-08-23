
rawType = 'tif';
neuronProbabilityType = 'png';
membraneProbabilityType = 'tiff';
mitoProbabilityType = 'png';

inputPath = '/home/thanuja/projects/toyData/set8/';
outputPath = '/home/thanuja/Dropbox/RESULTS/contourdetection/batch20140823/';
outputPathPNG = '/home/thanuja/Dropbox/RESULTS/contourdetection/batch20140823_png/';
% read all images in the raw images file path
rawImagePath = fullfile(inputPath,'raw');
allRawFiles = dir(fullfile(rawImagePath,'*.tif'));

% for each file
numFiles = length(allRawFiles);
for i=1:numFiles
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
end