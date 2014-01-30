function imageGrid = getImageGrid(imageIn,ofr,verbose)

% Inputs: 
%   imageIn: input image
%   ofr: oriented edge filter response

% Output:
%   ws_grid_equivalent = imageGrid
%   nodeEdges
%   edges2nodes
%   edges2pixels
%   nodeInds

% Parameters
gridResolution = 6;     % pixels

[sizeR,sizeC] = size(imageIn);

imageGrid = zeros(sizeR,sizeC);

%% node layout - square grid
% define nodes
[nodePix,numGridsX,numGridsY] = getGridNodeLayout_sq(sizeR,sizeC,gridResolution);
% nodePix: matrix containing the nodesPixInds in meshgrid format.
%% edge layout (node neighborhood) - square grid
% define edges2nodes
[adjacencyMat,nodeIndsVect] = getNodeAdjacency(nodePix);
[edges2nodes,nodeEdges] = getEdges2nodes_grid(adjacencyMat);
edges2pixels = getEdges2pixels_grid(edges2nodes,nodeIndsVect,sizeR,sizeC);
ws_grid = getWSfromGrid(nodeIndsVect,edges2pixels,nodeEdges,sizeR,sizeC,...
                gridResolution,numGridsX,numGridsY);
% for each node define tt

% edges2pixels ?

% nodeEdges: edgeIDs connected to each nodeID

