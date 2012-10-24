% train MRF from Dictionary patches

%% parameters 
bb = 16;
imgSize = 256;
outfilenameH = 'associationsH_1.mat';
outfilenameV = 'associationsV_1.mat';
%% load Dictionary
pathToDict = '/home/thanuja/matlabprojects/MRF_dict/Dictionary_256x_400w_16bb_NN.mat';
load(pathToDict);       % loads the structure output (from KSVD)
Dictionary = output.D;  % copy the dictionary
clear output;           % output from KSVD is no more required

%% get sparse coefficients for image 1
pathToSparseCoefMat = '/home/thanuja/matlabprojects/MRF_dict/sparse_stem1_256x_400w_16bb_NN.mat';
load(pathToSparseCoefMat);

%% get horizontal associations matrix
AMat_H = getHorizontalAssociations(sparsecoeff);

%% get vertical associations matrix
rowSize = imgSize -bb + 1;      % num of patches per row in image
AMat_V = getVerticalAssociations(sparsecoeff,rowSize);

%% sample MRF