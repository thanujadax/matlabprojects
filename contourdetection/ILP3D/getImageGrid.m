function [ws_grid,edgeSetRegions,edges2pixels,edges2nodes,nodeEdges,...
            adjacencyMat,nodeIndsVect,edges2regions,boundaryEdgeIDs,...
            twoRegionEdges]...
                                        = getImageGrid(imageIn,gridResolution,verbose)

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
% gridResolution = 6;     % pixels

[sizeR,sizeC] = size(imageIn);

%% node layout - square grid
% define nodes
[nodePix,numGridsX,numGridsY] = getGridNodeLayout_sq(sizeR,sizeC,gridResolution);
% nodePix: matrix containing the nodesPixInds in meshgrid format.
%% edge layout (node neighborhood) - square grid

[adjacencyMat,nodeIndsVect] = getNodeAdjacency_sq(nodePix);

[edges2nodes,nodeEdges] = getEdges2nodes_grid(adjacencyMat);

edges2pixels = getEdges2pixels_grid(edges2nodes,nodeIndsVect,sizeR,sizeC,gridResolution);

[ws_grid,edgeSetRegions,edges2regions] = getWSfromGrid_sq(nodePix,nodeIndsVect,...
                edges2pixels,nodeEdges,sizeR,sizeC,numGridsX,numGridsY);
            
boundaryEdgeIDs = getBoundaryEdgeIDs(edges2regions);

numEdges = size(edges2nodes,1);
allEdgesSeq = 1:numEdges;
twoRegionEdges = setdiff(allEdgesSeq,boundaryEdgeIDs);

%% Project OFR on to the grid

%% Supplementary functions

function boundaryEdgeIDs = getBoundaryEdgeIDs(edges2regions)
% edges bordering region 0 (col2) are boundary edges
boundaryEdgeIDs = find(edges2regions(:,2)==0);