function [localHoughSpaces,patchLocations rho, theta] = getLocalHoughSpaces(imgIn,rhoResolution,...
                        theta,bb,slidingDist)
                    
% for each patch
% calculate the hough space
% store the result in a cell array where each cell is a local Hough space

% output:
% localHoughTransforms: cell array of dimension [numRowPatches numColPatches].
%   the (i,j)th cell corresponds to the (i,j)th patch 
%   of size bb with overlapping sliding distance given by slidingDist.
%   Each cell contains a structure s where:
%       s.H - Hough space for the patch
%       s.origin - (row,col) of the origin relative to the entire image


imgSize = size(imgIn);
numRowPatch = (imgSize(2) - bb)/slidingDist + 1; % num patches per row
numColPatch = (imgSize(1) - bb)/slidingDist + 1; % num patches per column

% initialize cell array
localHoughSpaces = cell(numColPatch,numRowPatch);
patchLocations = cell(numColPatch,numRowPatch);
% initialize struct
%s = struct('H',houghSpace,'origin',[patchOriginY patchOriginX]);

% TODO: parallelize by storing the starting points in a different data
% structure
for i=1:numColPatch
    for j=1:numRowPatch
        patchOriginX = (mod(j,numRowPatch)-1)*slidingDist + 1;
        patchOriginY = (mod(i,numColPatch)-1)*slidingDist + 1;
        
        if(patchOriginY+bb-1<=imgSize(1) && patchOriginX+bb-1<=imgSize(2))
            localPatch = imgIn(patchOriginY:patchOriginY+bb-1,...
                                        patchOriginX:patchOriginX+bb-1);
            [houghSpace,theta,rho] = houghFixedLength(localPatch,rhoResolution,theta);
            
            localHoughSpaces{i,j} = houghSpace;
            patchLocations{i,j} = [patchOriginY patchOriginX];
        else
            localHoughSpaces{i,j} = 0;
        end
    end 
end

