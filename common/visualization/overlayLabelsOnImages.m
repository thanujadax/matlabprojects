function h = overlayLabelsOnImages(filePaths)

% For all images in directory

% Inputs:
%   backgroundImages - e.g. raw EM image. should be properly normalized
%   since we use imshow to display the images.
%   foregroundImages - e.g. segmentation to visualize overlaid on EM image


filePaths.outputPath = '/home/thanuja/projects/data/myelin/ssSEM/s909/overlaid';
filePaths.backgroundImagesPath = '/home/thanuja/projects/data/myelin/ssSEM/s909/raw';
filePaths.foregroundImagesPath = '/home/thanuja/projects/data/myelin/ssSEM/s909/output';
filePaths.fileType = 'png';

inputForegroundFilesString = strcat('*.',filePaths.fileType);
inputForegroundFilesString = fullfile...
    (filePaths.foregroundImagesPath,inputForegroundFilesString);
inputForegroundDirectory = dir(inputForegroundFilesString);

inputBackgroundFilesString = strcat('*.',filePaths.fileType);
inputBackgroundFilesString = fullfile...
    (filePaths.backgroundImagesPath,inputBackgroundFilesString);
inputBackgroundDirectory = dir(inputBackgroundFilesString);

numForegroundImages = length(inputForegroundDirectory);
numBackgroundImages = length(inputBackgroundDirectory);

str1 = sprintf('Num background images found: %d', numBackgroundImages);
disp(str1);
str1 = sprintf('Num foreground images found: %d', numForegroundImages);
disp(str1);

if(numBackgroundImages==numForegroundImages)
    numImagesToProcess = numBackgroundImages;
else
    numImagesToProcess = min(numBackgroundImages,numForegroundImages);
    disp('WARNING: numBackgroundImages ~= numForegroundImages!')
end

for i = 1:numImagesToProcess
    str1 = sprintf('Processing image pair %d out of %d',i,numImagesToProcess);
    disp(str1)
    
    backgroundImageFileName = fullfile...
        (filePaths.backgroundImagesPath,inputBackgroundDirectory(i).name);
    backgroundImage = imread(backgroundImageFileName);
    
    % kind of normalizing so that it's properly visualized using imshow
    backgroundImage = invertImage(backgroundImage);
    backgroundImage = backgroundImage./255;

    foregroundImageFileName = fullfile...
        (filePaths.foregroundImagesPath,inputForegroundDirectory(i).name);
    foregroundImage = imread(foregroundImageFileName);    

    figure;
    imshow(backgroundImage, 'InitialMag', 'fit')
    % Make a truecolor all-green image.
    green = cat(3, zeros(size(backgroundImage)),ones(size(backgroundImage)), zeros(size(backgroundImage)));
    hold 
    h = imshow(green); 

    hold off

    % Use our influence map as the 
    % AlphaData for the solid green image.
    foregroundImage = double(foregroundImage);
    foregroundImage(foregroundImage>0) = 0.1;
    set(h, 'AlphaData', foregroundImage)

    set(gca,'position',[0 0 1 1],'units','normalized'); % to remove borders of fig

    % save figure
    saveFileName = strcat('overlaid_',inputForegroundDirectory(i).name);
    saveFileName = fullfile(filePaths.outputPath,saveFileName);
    str1 = strcat('saving overlay image ',saveFileName);
    disp(str1);
    print(saveFileName,'-dpng');
end