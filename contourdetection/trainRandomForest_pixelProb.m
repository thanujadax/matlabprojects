function forest_pixelProb = trainRandomForest_pixelProb()

% Script to train Random Forest for pixel classification
%   1 - cell interior (bright)
%   0 - membrane (dark)

% parameters
pathForImages = '';
pathForLabels = '';

% read images
imgFiles = dir(pathForImages);
labelFiles = dir(pathForLabels);

% extract features
for i=1:length(imgFiles)
    imageFeatures = getFeatures();
     
end

% train forest