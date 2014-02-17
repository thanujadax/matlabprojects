function [nodePix,numEdgesX,numEdgesY] = getGridNodeLayout_sq(sizeR,sizeC,gridResolution)

% Inputs:
%   gridResolution - num pixels between 2 adjacenct nodes (edge length)

% Output:
%   nodePix: matrix containing the nodesPixInds in meshgrid format.


% define a square grid layout spanning the image and return the node
% pixels and the edges according to 4-neighborhood

numEdgesX = floor((sizeC-1)/(gridResolution));
marginPix_X = mod((sizeC-1),(gridResolution));
gridStartX = floor(marginPix_X/2);
gridStartX = max(1,gridStartX);

numEdgesY = floor((sizeR-1)/(gridResolution));
marginPix_Y = mod((sizeC-1),(gridResolution));
gridStartY = floor(marginPix_Y/2);
gridStartY = max(1,gridStartY);

% define node positions
x_pos = gridStartX:gridResolution:sizeC;
y_pos = gridStartY:gridResolution:sizeR;

[nodeX, nodeY] = meshgrid(x_pos,y_pos);

nodePix = sub2ind([sizeR sizeC],nodeY,nodeX);

% numGrids = numEdgesX * numEdgesY;