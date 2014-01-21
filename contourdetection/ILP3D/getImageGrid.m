function imageGrid = getImageGrid(imageIn,ofr)

% Inputs: 
%   imageIn: input image
%   ofr: oriented edge filter response

% Output:
%   image

% Parameters
gridResolution = 3;     % pixels

[sizeR,sizeC] = size(imageIn);

imageGrid = zeros(sizeR,sizeC);

% define nodes
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

% define edges2nodes

% for each node define tt

% edges2pixels ?

% nodeEdges: edgeIDs connected to each nodeID

