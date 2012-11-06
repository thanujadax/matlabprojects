% take in patch label vector and generate an image of the supposedly
% overlapping image patches as a collection of non-overlapping image
% patches, in order of occurence relative to the actual reconstruction

function [patchMapImage, patchMapMat] = getPatchMap(labelVector,Dictionary,dimX,dimY,bb)

% Input:
% labelVector - contains the label (word index as per dictionary) for each
% patch of the image as a row vector
% Dictionary - contains the words as column vectors
% dimX - number of pixels in the reconstructed image in X direction
% dimY - number of pixels in the reconstructed image in Y direction
% bb - block size

% Outputs:
% patchMapImage - as image (matrix containing the pixel values)
% patchMapMat - as a matrix containing the labels (word id)

%% 

% generate patch matrix using label info
numPatchPerRow = dimX - bb + 1;
numPatchPerCol = dimY - bb + 1;
patchMapMat = reshape(labelVector,numPatchPerRow,numPatchPerCol);

% initialize patch map image
dimXpatchMap = numPatchPerRow * bb;
dimYpatchMap = numPatchPerCol * bb;
patchMapMat = zeros(dimXpatchMap,dimYpatchMap);

% check if the inputs make sense
if(size(labelVector,2)~=numPatchPerRow*numPatchPerCol)
    display('ERROR:number of labels do not match the number of output patches!');
end

% fill the patch map image with the words of the dictionary
patchIndex = 1;
for i = 1:numPatchPerRow
    startingPoint_i = (i-1)*bb + 1;
    for j = 1:numPatchPerCol
        startingPoint_j = (j-1)*bb + 1;
        patchMapMat(startingPoint_i:(startingPoint_i+bb-1),...
            startingPoint_j:(startingPoint_j+bb-1)) = ...
                reshape(Dictionary(:,labelVector(patchIndex)),bb,bb);
        patchIndex = patchIndex + 1;
        
    end
end
