function [localHoughSpaces,patchLocations, rho, theta,maxHoughPeak] = getLocalHoughSpaces(imgIn,rhoResolution,...
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

% imgIn = flipud(imgIn);  % flip upside down. therefore, flipping is disabled in houghFixedLength ~ hough2.
imgSize = size(imgIn);
cols = floor((imgSize(2) - slidingDist)/(bb-slidingDist)); % num patches per row
rows = floor((imgSize(1) - slidingDist)/(bb-slidingDist)); % num patches per column

% initialize cell array
localHoughSpaces = cell(rows,cols);
patchLocations = cell(rows,cols);
% initialize struct
%s = struct('H',houghSpace,'origin',[patchOriginY patchOriginX]);

maxHoughPeak = 0; % init

% TODO: parallelize by storing the starting points in a different data
% structure
for i=1:rows
    patchOriginY = (i-1)*(bb-slidingDist) + 1;
    if(patchOriginY+bb-1<=imgSize(1))
        for j=1:cols
            patchOriginX = (j-1)*(bb-slidingDist) + 1;

            if(patchOriginX+bb-1<=imgSize(2))
                localPatch = imgIn(patchOriginY:patchOriginY+bb-1,...
                                            patchOriginX:patchOriginX+bb-1);
                %[houghSpace,theta,rho] = houghFixedLength(localPatch,rhoResolution,theta);
                [houghSpace,theta,rho] = hough(localPatch,'RhoResolution',rhoResolution,'Theta',theta);

                localHoughSpaces{i,j} = houghSpace;
                maxVal = max(max(houghSpace));
                if(maxVal>maxHoughPeak)
                    maxHoughPeak = maxVal;
                end
                patchLocations{i,j} = [patchOriginY patchOriginX];
            else
                localHoughSpaces{i,j} = 0;
            end
        end
    end
end

