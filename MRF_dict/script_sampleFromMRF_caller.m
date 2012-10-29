% script to call sampleFromMRF() - Gibbs sampler


bb = 16;
imgx = 48;
imgy = 48;
rowSize = imgx - bb +1;
colSize = imgy - bb +1;
numPatches = rowSize*colSize;
%initLabels = zeros(1,rowSize*colSize);

initLabels = ceil(rand(1,numPatches) * size(Dictionary,2));

inputData = zeros(1,size(initLabels,2));

% Dictionary
load('Dictionary_wellTrained.mat');

% verticalAMat
load('VerticalAssociations');
% horizontalAMat
load('HorizontalAssociations.mat');

% smoothening parameter for unary potentials (rms pixel error)
sigma = 1.2;
                
coefMat = sampleFromMRF(initLabels,inputData,Dictionary,rowSize,...
                    colSize,verticalA,horizontalA,sigma);