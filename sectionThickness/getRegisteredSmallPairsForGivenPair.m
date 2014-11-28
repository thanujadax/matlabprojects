function cell_registeredSmallPairs = getRegisteredSmallPairsForGivenPair(image1,image2,...
        patchSizeX,patchSizeY)
    
% Inputs:
% image1 - image matrix
% image2 - image matrix
% patchSizeX - size of small patches to be registered
% patchSizeY - size of small patches to be registered

% cell_registeredSmallPairs - cell array. each cell contains a registered pair of
% image matrices each of which are from the adjacent images.

% extract patches from image1

% register each patch with image2 and extract corresponding patches from
% image2
