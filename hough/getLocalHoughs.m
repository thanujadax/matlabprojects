function [localizedHoughTransforms, rho, theta] = getLocalHoughs(imgIn,rhoResolution,...
                        theta,bb,slidingDist)
                    
% for each patch
% calculate the hough space
% store the result in a cell arrary where each cell is a local Hough space

imgSize = size(imgIn);
numRowPatch = (imgSize(2) - bb)/slidingDist + 1; % num patches per row
numColPatch = (imgSize(1) - bb)/slidingDist + 1; % num patches per column

% initialize cell array
localizedHoughTransforms = cell(numRowPatch,numColPatch);


% TODO: parallelize by storing the starting points in a different data
% structure
for i=1:numColPatch
    for j=1:numRowPatch
        patchOriginX = (mod(j,numRowPatch)-1)*slidingDist + 1;
        patchOriginY = (mod(i,numColPatch)-1)*slidingDist + 1;
        
        if(patchOriginY+bb-1<=imgSize(1) && patchOriginX+bb-1<=imgSize(2))
            localPatch = imgIn(patchOriginY:patchOriginY+bb-1,patchOriginX:patchOriginX+bb-1);
            [houghSpace,theta,rho] = houghFixedLength(theImage,rhoResolution,theta);
            s = struct('H',houghSpace,'origin',[patchOriginY patchOriginX])
            localizedHoughTransforms{i,j} = s;
        else
            localizedHoughTransforms{i,j} = 0;
        end
    end 
end

