% script to call sampleFromMRF() - Gibbs sampler


bb = 16;
imgx = 48;
imgy = 48;

initLabels

% Dictionary
load('Dictionary_wellTrained.mat');

rowSize = imgx - bb +1;
colSize = imgy - bb +1;

% verticalAMat
load('VerticalAssociations');
% horizontalAMat
load('HorizontalAssociations.mat');

% smoothening parameter for unary potentials (rms pixel error)
sigma = 1.2;
                
coefMat = sampleFromMRF(initLabels,Dictionary,rowSize,...
                    colSize,verticalAMat,horizontalAMat,sigma);