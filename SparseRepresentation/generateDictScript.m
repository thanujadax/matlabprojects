% script calling generateDictionary
% Thanuja

%% Parameters
bb = 16;
K = 400;
maxNumBlocksToTrainOn = 1000; % training
maxBlocksToConsider = 1000;   % to adjust the sliding distance for reconstruction
L = 1;                          % maximum number of atoms per signal
sigma = 40;                     % error. useful when L is not set.
imageIn = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
slidingDis = 1;
numIterOfKsvd = 10;
C = 2;                          % factor for error. not useful when L is set. (>0)
NN = 1;
reduceDC = 0;
numBPiterations = 20;

[Dictionary output] = generateDictionary(bb,K,maxNumBlocksToTrainOn,...
    maxBlocksToConsider,sigma,imageIn, slidingDis,numIterOfKsvd,C,NN,L,reduceDC,numBPiterations);

