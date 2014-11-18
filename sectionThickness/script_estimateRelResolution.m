% script to run thickness estimation

inputImageFileName = '/home/thanuja/Dropbox/data/fibsem/smallSet2/raw/2.png';
% inputImageFileName = '/home/thanuja/projects/inputData/trainingHalf/raw/00.png';

inputImageMat1 = double(imread(inputImageFileName));
inputImageMat1 = inputImageMat1./255;

inputImageFileName = '/home/thanuja/Dropbox/data/fibsem/smallSet2/raw/3.png';
% inputImageFileName = '/home/thanuja/projects/inputData/trainingHalf/raw/01.png';

inputImageMat2 = double(imread(inputImageFileName));
inputImageMat2 = inputImageMat2./255;

[numR,numC] = size(inputImageMat1);

Imat = zeros(numR,numC,2);
Imat(:,:,1) = inputImageMat1;
Imat(:,:,2) = inputImageMat2;

maxG = 15; % maximum sampling distance to evaluate

[g,m12] = estimateRelResolution_xcorr(Imat,maxG);