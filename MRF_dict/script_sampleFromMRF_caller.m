% script to call sampleFromMRF() - Gibbs sampler

numIterations = 1;
sigma = 0.5;
bb = 16;
imgx = 256;
imgy = 256;
rowSize = imgx - bb +1;
colSize = imgy - bb +1;
numPatches = rowSize*colSize;
%initLabels = zeros(1,rowSize*colSize);

lambda = 1;  % weighting for the unary potential

maskSize = 10;
sigmaLPF = 5;

initLabels = ceil(rand(1,numPatches) * size(Dictionary,2));

%% input data

%filename = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
%unseen image
filename = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
inputImage = getImageIntoMatrix(filename);
inputImage = gaussianFilter(inputImage,maskSize,sigmaLPF);
inputData = im2col(inputImage(1:imgx,1:imgy),[bb,bb],'sliding');
% Dictionary
load('Dictionary_wellTrained_25w_scaled.mat');
% verticalAMat
load('VerticalAssociations_25w.mat');
% horizontalAMat
load('HorizontalAssociations_25w.mat');

% smoothening parameter for unary potentials (rms pixel error)
%%
                
[coefMat,labelVector] = sampleFromMRF(initLabels,inputData,Dictionary,rowSize,...
                    verticalA,horizontalA,sigma,lambda);
                
for i = 1:numIterations-1
    [coefMat,labelVector] = sampleFromMRF(labelVector,inputData,Dictionary,rowSize,...
                    verticalA,horizontalA,sigma,lambda);
end

%% Combine the overlapping patches to get the sampled whole image
IOut = getImageFromCoeff(Dictionary,coefMat,[imgx imgy],bb);
figure(3)
imagesc(IOut);
colormap('gray');


