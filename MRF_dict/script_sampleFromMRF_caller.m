% script to call sampleFromMRF() - Gibbs sampler

numIterations = 2;
sigma = 1;
bb = 16;
imgx = 256;
imgy = 256;
rowSize = imgx - bb +1;
colSize = imgy - bb +1;
numPatches = rowSize*colSize;
%initLabels = zeros(1,rowSize*colSize);

initLabels = ceil(rand(1,numPatches) * size(Dictionary,2));

%% input data

%filename = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';

filename = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
inputImage = getImageIntoMatrix(filename);
inputData = im2col(inputImage(1:imgx,1:imgy),[bb,bb],'sliding');

%%
% Dictionary
load('Dictionary_wellTrained.mat');

% verticalAMat
load('VerticalAssociations');
% horizontalAMat
load('HorizontalAssociations.mat');

% smoothening parameter for unary potentials (rms pixel error)

                
[coefMat,labelVector] = sampleFromMRF(initLabels,inputData,Dictionary,rowSize,...
                    colSize,verticalA,horizontalA,sigma);
                
for i = 1:numIterations-1
    [coefMat,labelVector] = sampleFromMRF(labelVector,inputData,Dictionary,rowSize,...
                    colSize,verticalA,horizontalA,sigma);
end

%% Combine the overlapping patches to get the sampled whole image
IOut = getImageFromCoeff(Dictionary,coefMat,[imgx imgy],bb);
figure(3)
imagesc(IOut);
colormap('gray');


