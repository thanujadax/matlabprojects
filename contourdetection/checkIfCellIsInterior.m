function isInterior = checkIfCellIsInterior(internalPixels,edgePixels,...
                    edgeOrientation,sizeR,sizeC)
% checks if a cell is part of cell interior based on the input edge
% orientation, which bounds the cell in concern
% returns a value between +1 and -1. If it's +1, it's most likely cell
% interior

% get center pixel of edge
numEdgePixels = numel(edgePixels);
centerpos = floor(numEdgePixels/2);
centerPixInd = edgePixels(centerpos);
[centerPix.y,centerPix.x] = ind2sub([sizeR sizeC], centerPixInd);

% determine the angle cog of cell makes with (wrt) the center pixel
numInternalPixels = numel(internalPixels);
isInterior = zeros(numInternalPixels,1);
for i=1:numInternalPixels
    theta = atan2d((internalPixels(i).y-centerPix.y),(internalPixels(i).x-centerPix.x));
    % convert theta into the same scale as edgeOrientation angles
    theta_scaled = convertThetaIntoOrientationScale(theta);

    % based on the angle and the edge orientation determine if the cell is interior or not
    isInterior(i) = sind(edgeOrientation - theta_scaled);

end

meanInteriorScore = mean(isInterior);
if(meanInteriorScore<0)
    isInterior = -1;
else
    isInterior = 1;
end
    

end

function theta2 = convertThetaIntoOrientationScale(theta)
if(theta>0 && theta <180)
    theta2 = 360 - theta;
elseif(theta<0 && theta>(-180))
    theta2 = -1 * theta;
else
    theta2 = theta;
end
    
end