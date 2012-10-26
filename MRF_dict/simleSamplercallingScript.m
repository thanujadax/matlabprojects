% script to run the sampleFromMRF function
imgDim = [256 256];
bb = 16;
load('Dictionary_wellTrained.mat');     % Dictionary
load('HorizontalAssociations.mat');     % horizontalA
load('VerticalAssociations');           % verticalA


[IOut,coefMat] = sampleFromMRF_simple(imgDim,bb,horizontalA,verticalA,Dictionary);