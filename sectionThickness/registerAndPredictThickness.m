function [t_median,t_mean,t_var] = registerAndPredictThickness(...
    inputImagePath,imgFileType,patchSizeX,patchSizeY)

% Inputs:

% Outputs:

% It could be difficult for two large consecutive serial section images
% to be globally registered properly. However, smaller regions of both
% images could be registered more effectively. Since the distance measure
% between images is based on correlation of coefficient (or related), such
% problems will introduce an unnecessary distance between two adjacent
% images. Therefore, for a given pair of images, we first register smaller
% sub-pairs of images and calculate the distance between each such pair and
% then take the median (to be robust against outliers)

allImageFiles = dir(fullfile(inputImagePath,strcat('*.',imgFileType)));
numImg = length(allImageFiles);
% for each pair of images
t_median = zeros(1,numImg-1);
t_mean = zeros(1,numImg-1);
t_var = zeros(1,numImg-1);
for i=1:numImg-1
    % extract corresponding smaller pieces and register
    image1 = fullfile(inputImageDir,allImageFiles(i).name);
    image2 = fullfile(inputImageDir,allImageFiles(i+1).name);
    cell_registeredSmallPairs = getRegisteredSmallPairs(image1,image2,...
        patchSizeX,patchSizeY);
    % calculate the distance between each pair
    estimatesForPair = getThicknessForRegisteredPairs(cell_registeredSmallPairs);
    [t_median(i),t_mean(i),t_var(i)] = thicknessStatsFromSmallPairs(estimatesForPair);

end
% return median mean and var

