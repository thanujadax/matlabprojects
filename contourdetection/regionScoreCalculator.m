function regionScores = regionScoreCalculator(forest,imIn)

%% perform cell interior classification - probability map

% create feature matrix for the input image
% Using membrane_RF/membraneDetection/membraneFeatures.m to extract
% features with rotational invariance
% The output is an n-by-n-by-14 3D matrix

% param for membrane-feature-extractor 
cs = 29;
ms = 3;
csHist = cs;

featureMat  = membraneFeatures(imIn, cs, ms, csHist);


% predict pixel probabilities using already trained Random Forest

% y_hat has the binary yes/no decision of the classifier
% votes is the votes of the trees (more like probability map)
[y_hat,pixelProbabilities] = classRF_predict(double(featureMat), forest);

%% calculate cell scores (mean, median, number of cells)

