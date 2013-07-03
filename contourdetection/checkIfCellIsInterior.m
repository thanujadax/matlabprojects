function isInterior = checkIfCellIsInterior(internalPixels,edgePixels,...
                    edgeOrientation,sizeR,sizeC)
% checks if a cell is part of cell interior based on the input edge
% orientation, which bounds the cell in concern
% returns a value between +1 and -1. If it's +1, it's most likely cell
% interior

% Inputs:
%   internalPixels: matrix of x y coordinates. col1:x, col2:y

% get center pixel of edge
numEdgePixels = numel(edgePixels);
centerpos = floor(numEdgePixels/2);
if(centerpos==0);
    centerpos=1;
end
centerPixInd = edgePixels(centerpos);
[centerPix.y,centerPix.x] = ind2sub([sizeR sizeC], centerPixInd);

% determine the angle cog of cell makes with (wrt) the center pixel
numInternalPixels = size(internalPixels,1);
isInterior = zeros(numInternalPixels,1);

theta = atan2d((internalPixels(:,2)-centerPix.y),(internalPixels(:,1)-centerPix.x));
% convert theta into the same scale as edgeOrientation angles
theta_scaled = convertThetaIntoOrientationScale(theta);

% based on the angle and the edge orientation determine if the cell is interior or not
isInterior = sind(edgeOrientation - theta_scaled);


maxInteriorScore = max(isInterior);

if(maxInteriorScore>0)
    isInterior = 1;
else
    isInterior = - 1;
end
    
end

function theta2 = convertThetaIntoOrientationScale(theta)
numTheta = numel(theta);
theta2 = theta; % initialization
for i=1:numTheta
    if(theta(i)>0 && theta(i) <180)
        theta2(i) = 360 - theta(i);
    elseif(theta(i)<0 && theta(i)>(-180))
        theta2(i) = -1 * theta(i);
    end
end    
end