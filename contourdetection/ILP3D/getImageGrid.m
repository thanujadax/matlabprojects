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
gridResolution = 4;     % pixels

[sizeR,sizeC] = size(imageIn);

imageGrid = zeros(sizeR,sizeC);

%% node layout - square grid
% define nodes
nodePix = getGridNodeLayout_sq(sizeR,sizeC,gridResolution);

%% edge layout (node neighborhood) - square grid
% define edges2nodes
adjacencyMat = getNodeAdjacency(nodePix);
edges2nodes = getEdges2nodes_grid(adjacencyMat);
edges2pixels = getEdges2pixels_grid();
% for each node define tt

% edges2pixels ?

% nodeEdges: edgeIDs connected to each nodeID

