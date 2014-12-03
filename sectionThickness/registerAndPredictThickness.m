function [t_median,t_mean,t_std, thicknessCurve,tCurve_std] ...
    = registerAndPredictThickness(...
    inputImagePath,imgFileType,patchSizeX,patchSizeY,maxNumPatches,overlap,...
    maxThicknessPix)

% Inputs:

% maxNumPatches - maximum number of small patches per image to be
% considered

% Outputs:

% It could be difficult for two large consecutive serial section images
% to be globally registered properly. However, smaller regions of both
% images could be registered more effectively. Since the distance measure
% between images is based on correlation of coefficient (or related), such
% problems will introduce an unnecessary distance between two adjacent
% images. Therefore, for a given pair of images, we first register smaller
% sub-pairs of images and calculate the distance between each such pair and
% then take the median (to be robust against outliers)

% NB: If the image patches are  too small it will negatively affect the
% thickness estimate. If the image patches are too large, the patch registration
% process will not be effective.

allImageFiles = dir(fullfile(inputImagePath,strcat('*.',imgFileType)));
numImg = length(allImageFiles);
% for each pair of images
t_median = zeros(1,numImg-1);
t_mean = zeros(1,numImg-1);
t_std = zeros(1,numImg-1);


[thicknessCurve, tCurve_std] = getThicknessCurve(inputImagePath,imgFileType,maxThicknessPix);

% plot estimator curve
figure();
shadedErrorBar((1:maxThicknessPix),thicknessCurve,tCurve_std,'g');

for i=1:numImg-1
    % extract corresponding smaller pieces and register
    image1 = fullfile(inputImageDir,allImageFiles(i).name);
    image2 = fullfile(inputImageDir,allImageFiles(i+1).name);
    [image1_patches,image2_patches] = getRegisteredSmallPairs(image1,image2,...
        patchSizeX,patchSizeY,maxNumPatches,overlap);
    % calculate the distance between each pair
    [t_median(i),t_mean(i),t_std(i)] = getThicknessForRegisteredPairs...
            (image1_patches,image2_patches,thicknessCurve,maxThicknessPix);
    

end
% plot estimated thickness
figure();
shadedErrorBar((1:maxThicknessPix),t_median,t_std,'b');
