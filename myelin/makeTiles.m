function makeTiles(imageRoot,outputPath,params)

% Inputs:
% imageRoot - root folder of subdirectories each containing large tif images
% outputPath - output tiles are written here with the same subdirectory
% hierarchy as the image root
% params.xBorder = 25; % px
% params.yBorder = 25; % px
% params.xTiles = 40;
% params.yTiles = 40;

%% Test inputs
imageRoot = '/home/thanuja/projects/data/ssSEM_dataset/tiles/80';
outputPath = '/home/thanuja/projects/data/ssSEM_dataset/tinyTiles/80';
params.xBorder = 25; % px
params.yBorder = 25; % px
params.xTiles = 20;
params.yTiles = 20;
params.imageType = 'tif';

%% Method description
% get all subdirectories in imageRoot
% in each subdirectory read the large tif images
% from it's size, calculate the grid points for all the regions of interest
% extract tiles adding the borders to the regions of interest
% save the images in outputPath with under the same subdir name as input
% image
% info file - roi_x and roi_y lengths in pixels for each image

% get the list of directories
[sampleDirectories,~] = subdir(imageRoot);

% tiff images are stored in each directory
parfor i=1:length(sampleDirectories)
    
    sampleSubDirName = sampleDirectories{i};
    % read all image stacks in this sample
    imageFileString = strcat('*.',params.imageType);
    imageFileString = fullfile(sampleSubDirName,imageFileString);
    imageSubDir = dir(imageFileString);
    
    str1 = sprintf('Making tiles from image %s',sampleSubDirName);
    disp(str1)
    % process each image stack in the sample
    for j=1:length(imageSubDir)
        inputImageFileName = fullfile...
            (sampleSubDirName,imageSubDir(j).name);
    tokenizedSubDirName = strsplit(sampleSubDirName,filesep);
    tokenizedSubDirName = tokenizedSubDirName{end};
    % make subdir to save if doesn't exist
    checkAndCreateSubDir(outputPath,tokenizedSubDirName);
    outputSavePath = fullfile(outputPath,tokenizedSubDirName);  
        
    makeTilesFromImage(inputImageFileName,params,outputSavePath)
    % writes output to output path as txt file. Col vector.
    end
    
    % TODO: do we use the same curve for ths same sample?
    % makes sense to do so. Make aggregate curve?
end