function regionScores = regionScoreCalculator(forest,imIn,setOfCells,edges2pixels,...
    nodeInds,edges2nodes,K,wsIndsForCells,ws,showImg)

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

[sizeR,sizeC] = size(imIn);

% TODO: modify feature extraction to use OFR
% TODO: check if feature matrix for this image is already calculated
fm  = membraneFeatures(imIn, maskSize, halfwidthOfStrElement, csHist);

% reshape feature matrix so that each column stands for a feature and each
% row corresponds to a pixel (object to be classified)
fm = reshape(fm,size(fm,1)*size(fm,2),size(fm,3));
fm(isnan(fm))=0;

% predict pixel probabilities using already trained Random Forest

% y_hat has the binary yes/no decision of the classifier
% votes is the votes of the trees (more like probability map)
[y_hat,v] = classRF_predict(double(fm), forest);

% visualize prediction
pixelProbabilities = v(:,2);

pixelProbabilities = reshape(pixelProbabilities,[sizeR sizeC]);
pixelProbabilities = double(pixelProbabilities)/max(pixelProbabilities(:));
% figure; imshow(makeColorOverlay(pixelProbabilities,imIn));
% title('predicted pixel probabilities overlay')
%% calculate region scores (mean, median, number of cells)
% if region score is close to 1 it's probably cell interior
% if close to zero, it's probably membrane/(mitochondria?)

% for each cell, get its interior pixels
% get the mean of the pixelProbabilities for these pixels
% TODO: how to use the size of the region?

% have to invert the probabilities. here, 1 correspond to dark pixels.
% change it around
pixelProbabilities = (pixelProbabilities-1).*(-1);
if(showImg)
    figure;imshow(pixelProbabilities);
    title('pixel probability map')
end
regionScores = getCellPriors_probability(pixelProbabilities,setOfCells,...
    sizeR,sizeC,wsIndsForCells,ws,showImg);


