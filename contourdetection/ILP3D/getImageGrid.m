function [ws_grid,edgeSetRegions,edges2pixels,edges2nodes,nodeEdges,...
            adjacencyMat,nodeIndsVect]...
                                        = getImageGrid(imageIn,verbose)

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

%% node layout - square grid
% define nodes
[nodePix,numGridsX,numGridsY] = getGridNodeLayout_sq(sizeR,sizeC,gridResolution);
% nodePix: matrix containing the nodesPixInds in meshgrid format.
%% edge layout (node neighborhood) - square grid
% define edges2nodes
[adjacencyMat,nodeIndsVect] = getNodeAdjacency(nodePix);
[edges2nodes,nodeEdges] = getEdges2nodes_grid(adjacencyMat);
edges2pixels = getEdges2pixels_grid(edges2nodes,nodeIndsVect,sizeR,sizeC);
[ws_grid, edgeSetRegions] = getWSfromGrid(nodeIndsVect,edges2pixels,nodeEdges,...
                sizeR,sizeC,numGridsX,numGridsY);

%% project OFR on to the grid