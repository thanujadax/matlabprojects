% function Dictionary = generateDictionary(bb,RR,K,maxNumBlocksToTrainOn,sigma,imageIn)
% Generates a dictionary using an input image
% removes Gaussian noise from the image
% uses KSVD Matlab Toolbox
%
% Thanuja 05.09.2012

function [Dictionary output] = generateDictionary(bb,RR,K,maxNumBlocksToTrainOn,...
    maxBlocksToConsider,sigma,imageIn, slidingDis,numIterOfKsvd,C,NN)

% NN - nonnegative flag

[IMin,pp]=imread(imageIn);

% fixing input image format
IMin=im2double(IMin);
if (length(size(IMin))>2)
    IMin = rgb2gray(IMin);
end
if (max(IMin(:))<2)
    IMin = IMin*255;
end

%%  KSVD routine to train dictionary
% Learn dictionary
reduceDC = 1;
[NN1,NN2] = size(IMin);
waitBarOn = 1;
displayFlag = 1;

% train a dictionary on blocks from the noisy image
if(prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn)
    randPermutation =  randperm(prod([NN1,NN2]-bb+1));
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);

    blkMatrix = zeros(bb^2,maxNumBlocksToTrainOn);
    for i = 1:maxNumBlocksToTrainOn
        [row,col] = ind2sub(size(IMin)-bb+1,selectedBlocks(i));
        currBlock = IMin(row:row+bb-1,col:col+bb-1);
        blkMatrix(:,i) = currBlock(:);
    end
else
    blkMatrix = im2col(IMin,[bb,bb],'sliding');
end

param.K = K;
param.numIteration = numIterOfKsvd ;

%param.errorFlag = 1; 
% decompose signals until a certain error is reached. do not use fix number of coefficients.

param.errorFlag = 0; % describe each signal by exactly one atom
param.L = 1;

param.errorGoal = sigma*C;
param.preserveDCAtom = 0;  % if 1, tells KSVD to keep the first word of Dict constant

Pn=ceil(sqrt(K));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

param.initialDictionary = DCT(:,1:param.K );
param.InitializationMethod =  'GivenMatrix';

if (reduceDC)
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end

if (waitBarOn)
    counterForWaitBar = param.numIteration+1;
    h = waitbar(0,'Denoising In Process ...');
    param.waitBarHandle = h;
    param.counterForWaitBar = counterForWaitBar;
end


param.displayProgress = displayFlag;
if(NN)
    [Dictionary,output] = KSVD_NN(blkMatrix,param);
else
    [Dictionary,output] = KSVD(blkMatrix,param);
end

