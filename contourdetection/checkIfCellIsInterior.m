function isInterior = checkIfCellIsInterior(cog,edgePixels,edgeOrientation,sizeR,sizeC)
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
theta = tan2d((cog.y-centerpixel.y),(cog.x-centerpixel.x));
% convert theta into the same scale as edgeOrientation angles
theta_scaled = convertThetaIntoOrientationScale(theta);

% based on the angle and the edge orientation determine if the cell is interior or not
isInterior = sind(edgeOrientation - theta_scaled);

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