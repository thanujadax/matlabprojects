% Algorithm that denoises unseen images using a denoised dictionary, for EM
% images
% uses the KSVD MATLAB Toolbox
% 
% Thanuja 05.09.2012

%% Parameters
LearnNewDictionaryN = 0;
LearnNewDictionary0 = 0;
bb = 16; % block size
RR = 4; % redundancy factor
%K = RR*bb^2; % number of atoms in the dictionary
K = 400;
maxNumBlocksToTrainOn = 1000; %  - the maximal number of blocks
%                       to train on. The default value for this parameter is
%                       65000. However, it might not be enough for very large
%                       images
maxBlocksToConsider = 60000; % - maximal number of blocks that
%                       can be processed. This number is dependent on the memory
%                       capabilities of the machine, and performances�
%                       considerations. If the number of available blocks in the
%                       image is larger than 'maxBlocksToConsider', the sliding
%                       distance between the blocks increases. The default value
numWords = 1;       % every patch should be represented by exactly one atom.
                    % if set to zero, it keeps the number of atoms open and
                    % keeps within the error goal.
sigma = 80;
pathForImages = '/home/thanuja/matlabprojects/data/mitoData/';
imageName = 'stem1_256by256.png'; % for initial dictionary

newImagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256_2.png'; % unseen image
IMin_0 = imread(newImagePath);
% preprocessing (to remove the dark edges at the bottom and to the right)
removedEdgeSize = 0;
if (1)
    imSize = size(IMin_0);
    resizeDim = [imSize(1)-removedEdgeSize imSize(2)-removedEdgeSize];
    IMin = IMin_0(1:resizeDim(1), 1:resizeDim(2), 1);
end
% fixing input image format
IMin=im2double(IMin);
if (length(size(IMin))>2)
    IMin = rgb2gray(IMin);
end
if (max(IMin(:))<2)
    IMin = IMin*255;
end

imageIn = strcat([pathForImages,imageName]);
C = 1.15;  % error factor: to control the weight on sigma for denoising the noisy image to learn the dictionary
slidingDis = 1; % the gap between two blocks considered in the sliding window
numIterOfKsvd = 5; % for training the dictionary - the number of KSVD iterations processed
%                       blocks from the noisy image. 
waitBarOn = 1;
reduceDC = 1;

% TEST: loads Dictionary (from noisy image)
load('/home/thanuja/matlabprojects/MRF_dict/Dictionary_256x_400w_16bb_NN.mat');
Dictionary = output.D;  % copy the dictionary from 'output' struct loaded
clear output;           % output from KSVD is no more required 
%% Learn the dictionary Dn from a noisy image In (if not already provided)
if(LearnNewDictionaryN)
[Dictionary_n, output_n] = generateDictionary(bb,RR,K,maxNumBlocksToTrainOn,...
    maxBlocksToConsider,sigma,imageIn, slidingDis,numIterOfKsvd,C);
end

%% Learn the dictionary D0 corresponding to the denoised version of image In
if(LearnNewDictionary0)
% Denoise Dn to get D0

end

%% Reconstruct the unseen image
[IOut,sparsecoeff,vecOfMeans] = OMPDenoisedImage(IMin,Dictionary,bb,...
    maxBlocksToConsider,sigma,C,slidingDis,waitBarOn,reduceDC,numWords);
imshow(IOut,[]);