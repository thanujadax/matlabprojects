function [image1_patches,image2_patches] = getRegisteredSmallPairsForGivenPair...
    (image1,image2,patchSizeX,patchSizeY,maxNumPatches,overlap)
    
% Inputs:
% image1 - image matrix
% image2 - image matrix
% patchSizeX - size of small patches to be registered
% patchSizeY - size of small patches to be registered

% cell_registeredSmallPairs - cell array. each cell contains a registered pair of
% image matrices each of which are from the adjacent images.

% extract patches from image1
image1_patches = extractPatchesFromImage(image1,patchSizeX,patchSizeY,...
                    maxNumPatches,overlap);

% register each patch with image2 and extract corresponding patches from
% image2
image2_patches = extractRegisteredPatchesForGivenPatches(...
                                    image2,image1_patches);
