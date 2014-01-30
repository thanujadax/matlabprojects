function [ws_grid,edgeSetRegions] = getWSfromGrid_sq(nodeIndsVect,edges2pixels,nodeEdges,sizeR,sizeC,...
                    gridResolution,numGridsX,numGridsY)
                
% Const
NUM_SIDES = 4; % number of sides per  grid cell - 4 for square shaped grid cells

% generates a ws like output from the square grid representation
% all edges and nodes will be labeled zero where as all other 'regions'
% will be given a unique ID, including the image border (=1)

% init
ws_grid = ones(sizeR,sizeC);
% border already assinged label '1'

% mark all edges and nodes with label '0'
ws_grid(nodeIndsVect) = 0;
% remove 1st col (edgeIDs)
edges2pixels(:,1) = [];
% get all edge pixel inds (nonzero)
allEdgePixelInds = edges2pixels(edges2pixels>0);
ws_grid(allEdgePixelInds) = 0;

% assign unique region IDs (starting from 2)
% (using the fact that the grid cells are squares)

numGridRegions = numGridsX * numGridsY;
% init
edgeSetRegions = zeros(numGridRegions,NUM_SIDES);
nodes4_tmp = zeros(4,1);
edgeSet_region_i = zeros(4,1);
k = 0;
for i=1:(numGridsY-1)
    for j=1:(numGridsX-1)
    % get the 4 nodes for grid_i
    k = k+1;
    nodes4_tmp(1) = nodePix(i,j);
    nodes4_tmp(2) = nodePix(i,(j+1));
    nodes4_tmp(3) = nodePix((i+1),(j+1));
    nodes4_tmp(4) = nodePix((i+1),j);
    % get the 4 edges using the node info -> edgeSetRegions
    nodeEdgeSets_tmp1 = nodeEdges(nodes4_tmp(1),:);
    nodeEdgeSets_tmp1 = nodeEdgeSets_tmp1(nodeEdgeSets_tmp1>0);
    nodeEdgeSets_tmp2 = nodeEdges(nodes4_tmp(2),:);
    nodeEdgeSets_tmp2 = nodeEdgeSets_tmp2(nodeEdgeSets_tmp2>0);
    nodeEdgeSets_tmp3 = nodeEdges(nodes4_tmp(3),:);
    nodeEdgeSets_tmp3 = nodeEdgeSets_tmp3(nodeEdgeSets_tmp3>0);
    nodeEdgeSets_tmp4 = nodeEdges(nodes4_tmp(4),:);
    nodeEdgeSets_tmp4 = nodeEdgeSets_tmp4(nodeEdgeSets_tmp4>0);
    % get region pixels -> ws_grid
    edgeSet_region_i(1) = intersect(nodeEdgeSets_tmp1,nodeEdgeSets_tmp2);
    edgeSet_region_i(2) = intersect(nodeEdgeSets_tmp2,nodeEdgeSets_tmp3);
    edgeSet_region_i(3) = intersect(nodeEdgeSets_tmp3,nodeEdgeSets_tmp4);
    edgeSet_region_i(4) = intersect(nodeEdgeSets_tmp4,nodeEdgeSets_tmp1);
    
    edgeSet_region_i = edgeSet_region_i(edgeSet_region_i>0);
    k = k+1;
    edgeSetRegions(k,1:4) = edgeSet_region_i;    
    
    % regionIDs (=k)
    % get the pixels for the grid region
    [rNodes,cNodes] = ind2sub([sizeR sizeC],nodes4_tmp);
    rStart = rNodes(1) + 1;
    rStop = rNodes(3) - 1;
    cStart = cNodes(1) + 1;
    cStop = cNodes(2) -1;
    
    % assign unique regionID to regionPixels_k (k+1). ID=1 is the border
    ws_grid(rStart:rStop,cStart:cStop) = k+1;
    
    end
end


