function nodePix = getGridNodeLayout_sq(sizeR,sizeC,gridResolution)

% Output:
%   nodePix: matrix containing the nodesPixInds in meshgrid format.


% define a square grid layout spanning the image and return the node
% pixels and the edges according to 4-neighborhood

numEdgesX = floor(sizeC/gridResolution);
marginPix_X = mod(sizeC,gridResolution);
gridStartX = floor(marginPix_X/2);

numEdgesY = floor(sizeR/gridResolution);
marginPix_Y = mod(sizeC,gridResolution);
gridStartY = floor(marginPix_Y/2);

% define node positions
x_pos = gridStartX:gridResolution:sizeC;
y_pos = gridStartY:gridResolution:sizeR;

[nodeX, nodeY] = meshgrid(x_pos,y_pos);

nodePix = sub2ind([sizeR sizeC],nodeY,nodeX);