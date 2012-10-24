% script calling generateDictionary
% Thanuja

%% Parameters
bb = 16;
K = 400;
maxNumBlocksToTrainOn = 65000;  % training
maxBlocksToConsider = 100000;   % to adjust the sliding distance for reconstruction
L = 1;                          % maximum number of atoms per signal
sigma = 40;                     % error. useful when L is not set.
imageIn = '/home/thanuja/matlabprojects/data/mitoData/stem1.png';
slidingDis = 1;
numIterOfKsvd = 10;
C = 2;                          % factor for error. not useful when L is set. (>0)
NN = 1;

[Dictionary output] = generateDictionary(bb,K,maxNumBlocksToTrainOn,...
    maxBlocksToConsider,sigma,imageIn, slidingDis,numIterOfKsvd,C,NN)

