% train MRF from Dictionary patches

% load Dictionary
pathToDict = 'Dictionary_256x_400w_16bb_NN.mat';
load(pathToDict);       % loads the structure output (from KSVD)
Dictionary = output.D;  % copy the dictionary
clear output;           % output from KSVD is no more required

% get sparse coefficients for image 1

% get horizontal associations matrix

% get vertical associations matrix