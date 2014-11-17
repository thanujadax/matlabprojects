% script to run thickness estimation sporring2014

inputImageFileName = '/home/thanuja/Dropbox/data/fibsem/smallSet2/raw/0.png';

inputImageMat1 = double(imread(inputImageFileName));
inputImageMat1 = inputImageMat1./255;

inputImageFileName = '/home/thanuja/Dropbox/data/fibsem/smallSet2/raw/1.png';

inputImageMat2 = double(imread(inputImageFileName));
inputImageMat2 = inputImageMat2./255;

[numR,numC] = size(inputImageMat1);

Imat = zeros(numR,numC,2);
Imat(:,:,1) = inputImageMat1;
Imat(:,:,2) = inputImageMat2;

maxG = 3; % maximum sampling distance to evaluate

g = estimateSamplingRatio(Imat,maxG)