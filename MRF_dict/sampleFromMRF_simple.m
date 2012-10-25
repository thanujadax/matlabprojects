%% function IOut = sampleFromMRF_simple(imgDim,bb,Hmat,Vmat,prior,Dictionary)

% generates a sample image based on the learned MRF. Generates the image
% sequentially from the top left corner to the bottom right corner.

% Inputs:
% imgDim - vector containing the size of the output image to be sampled [rows columns]
% bb - block size
% Hmat - horizontal association matrix learned (row,column -> left to right)
% Vmat - vertical association matrix learned (row,column -> top to bottom)
% prior - prior for label (word) usage
% Dictionary

% Output:
% Iout - sampled image

function IOut = sampleFromMRF_simple(imgDim,bb,Hmat,Vmat,prior,Dictionary) 

% Initialize random field
rows = imgDim(1) - bb + 1;              % number of overlapping image patches per column
cols = imgDim(2) - bb + 1;              % number of overlapping image patches per row

totPatches = rows*cols;                 % total number of patches to be generated

coefMat = sparse(size(Dictionary,2),totPatches);

% Pick a random word for the first image patch
wordInd = ceil(rand(1) * size(Dictionary,2)) ;
coefMat(wordInd,1) = 1;                 % assign the word to the first patch
for i = 2:totPatches
    % loop until all the patches are processed
    
    
end
