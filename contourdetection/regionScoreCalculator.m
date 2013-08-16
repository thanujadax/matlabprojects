function regionScores = regionScoreCalculator(forest,imIn,setOfCells,edges2pixels,...
    nodeInds,edges2nodes,K,sizeR,sizeC)

% Inputs:
%   imgIn - normalized image. 1 -> bright
%   K - positive scalar factor for the costs

%% perform cell interior classification - probability map

% create feature matrix for the input image
% Using membrane_RF/membraneDetection/membraneFeatures.m to extract
% features with rotational invariance
% The output is an n-by-n-by-14 3D matrix

% param for membrane-feature-extractor 
maskSize = 29;
halfwidthOfStrElement = 3; % actual width = halfwidth + 1
csHist = maskSize;

% TODO: modify feature extraction to use OFR
featureMat  = membraneFeatures(imIn, maskSize, halfwidthOfStrElement, csHist);


% predict pixel probabilities using already trained Random Forest

% y_hat has the binary yes/no decision of the classifier
% votes is the votes of the trees (more like probability map)
[y_hat,pixelProbabilities] = classRF_predict(double(featureMat), forest);

%% calculate region scores (mean, median, number of cells)
% if region score is close to 1 it's probably cell interior
% if close to zero, it's probably membrane/(mitochondria?)

% for each cell, get its interior pixels
% get the mean of the pixelProbabilities for these pixels
% TODO: how to use the size of the region?

regionScores = getCellPriors_probability(pixelProbabilities,setOfCells,edges2pixels,...
    nodeInds,edges2nodes,K,sizeR,sizeC);
