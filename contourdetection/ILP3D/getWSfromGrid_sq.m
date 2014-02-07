function [ws_grid,edgeSetRegions,edges2regions] = getWSfromGrid_sq...
                    (nodePix,nodeIndsVect,edges2pixels,...
                    nodeEdges,sizeR,sizeC,numGridsX,numGridsY)
                
% Outputs:
%   ws_grid: watershed like output for the image organized in a uniform
%   square grid
%   edges2regions: rowID = edgeID, col1= regionID, col2= regionID.
%   regionIDs start with 1 without considering the image border
%   edgeSetRegions: rowID = regionID, without considering the image border
% 
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

numEdges = size(edges2pixels,1);
edges2regions = zeros(numEdges,2);

% assign unique region IDs (starting from 2)
% (using the fact that the grid cells are squares)

numGridRegions = numGridsX * numGridsY;
% init
edgeSetRegions = zeros(numGridRegions,NUM_SIDES);
nodePixInds4_tmp = zeros(4,1);
nodeListInds4_tmp = zeros(4,1);

edgeSet_region_k = zeros(4,1);

% get rid of the first column of nodeEdges (it contains nodeInds)
nodeEdges(:,1) = [];
k = 0;
for i=1:numGridsY
    for j=1:numGridsX
    % get the 4 nodes for grid_i
    
    nodePixInds4_tmp(1) = nodePix(i,j);
    [~,nodeListInds4_tmp(1)] = intersect(nodeIndsVect,nodePixInds4_tmp(1));
    nodePixInds4_tmp(2) = nodePix(i,(j+1));
    [~,nodeListInds4_tmp(2)] = intersect(nodeIndsVect,nodePixInds4_tmp(2));
    nodePixInds4_tmp(3) = nodePix((i+1),(j+1));
    [~,nodeListInds4_tmp(3)] = intersect(nodeIndsVect,nodePixInds4_tmp(3));
    nodePixInds4_tmp(4) = nodePix((i+1),j);
    [~,nodeListInds4_tmp(4)] = intersect(nodeIndsVect,nodePixInds4_tmp(4));
    
    % get the 4 edges using the node info -> edgeSetRegions
    nodeEdgeSets_tmp1 = nodeEdges(nodeListInds4_tmp(1),:);
    nodeEdgeSets_tmp1 = nodeEdgeSets_tmp1(nodeEdgeSets_tmp1>0);
    nodeEdgeSets_tmp2 = nodeEdges(nodeListInds4_tmp(2),:);
    nodeEdgeSets_tmp2 = nodeEdgeSets_tmp2(nodeEdgeSets_tmp2>0);
    nodeEdgeSets_tmp3 = nodeEdges(nodeListInds4_tmp(3),:);
    nodeEdgeSets_tmp3 = nodeEdgeSets_tmp3(nodeEdgeSets_tmp3>0);
    nodeEdgeSets_tmp4 = nodeEdges(nodeListInds4_tmp(4),:);
    nodeEdgeSets_tmp4 = nodeEdgeSets_tmp4(nodeEdgeSets_tmp4>0);
    % get region pixels -> ws_grid
    edgeSet_region_k(1) = intersect(nodeEdgeSets_tmp1,nodeEdgeSets_tmp2);
    edgeSet_region_k(2) = intersect(nodeEdgeSets_tmp2,nodeEdgeSets_tmp3);
    edgeSet_region_k(3) = intersect(nodeEdgeSets_tmp3,nodeEdgeSets_tmp4);
    edgeSet_region_k(4) = intersect(nodeEdgeSets_tmp4,nodeEdgeSets_tmp1);
    
    edgeSet_region_k = edgeSet_region_k(edgeSet_region_k>0);
    k = k+1;
    edgeSetRegions(k,1:4) = edgeSet_region_k;    
    
    % update edges2regions (regionID = k)
    numEdgesForRegion_k = numel(edgeSet_region_k);
    for m=1:numEdgesForRegion_k
        edgeID_m = edgeSet_region_k(m);
        if(edges2regions(edgeID_m,1)==0)
            % 1st col is still empty. update 1st col
            edges2regions(edgeID_m,1) = k;
        else
            % 1st col is not empty. update 2nd col
            edges2regions(edgeID_m,2) = k;
        end
    end
    
    % regionIDs (=k)
    % get the pixels for the grid region
    [rNodes,cNodes] = ind2sub([sizeR sizeC],nodePixInds4_tmp);
    rStart = rNodes(1) + 1;
    rStop = rNodes(3) - 1;
    cStart = cNodes(1) + 1;
    cStop = cNodes(2) -1;
    
    % assign unique regionID to regionPixels_k (k+1). ID=1 is the border
    ws_grid(rStart:rStop,cStart:cStop) = k+1;
    
    end
end


