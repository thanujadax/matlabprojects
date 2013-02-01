% script calling generateDictionary
% Thanuja

%% Parameters
bb = 16;
K = 625; % number of words
maxNumBlocksToTrainOn = 100000; % training
maxBlocksToConsider = 100000;   % to adjust the sliding distance for reconstruction
L = 5;                          % maximum number of atoms per signal
sigma = 40;                     % error. useful when L is not set.
%imageIn = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
imageIn = '/home/thanuja/Dropbox/data/em_2013january/raw/00.tif';
slidingDis = 1;
numIterOfKsvd = 100;
C = 2;                          % factor for error. not useful when L is set. (>0)
NN = 1;
reduceDC = 0;
numBPiterations = 100;

[Dictionary output] = generateDictionary(bb,K,maxNumBlocksToTrainOn,...
    maxBlocksToConsider,sigma,imageIn, slidingDis,numIterOfKsvd,C,NN,L,reduceDC,numBPiterations);

