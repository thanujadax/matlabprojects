function c_edgesForRegions_cw = getOrderedRegionEdgeListIndsDir...
        (setOfRegions,edges2nodes_directional,nodeEdges,jAnglesAll_alpha)

% at each node take a right turn. get the edge. check if that edge is in
% the set of edges

numRegions = size(setOfRegions,1);
c_edgesForRegions_cw = cell(numRegions,1);

for i=1:numRegions
    edgesOfRegion_i = setOfRegions(i,:);
    edgesOfRegion_i = edgesOfRegion_i(edgesOfRegion_i>0);
    % arrange the edges in clockwise order of node traversal
    c_edgesForRegions_cw{i} = getCWorderedEdgesForRegion();
    
end


function cwOrderedEdgeListInds = getCWorderedEdgesForRegion...
                (edgeListInds_region,edges2nodes_directional)
% Pick one edge (1st in the list)
edgeListInd_1 = edgeListInds_region(1);
% Get the nodes at each end of the edge. At each node get the next edge as
% if to complete a clockwise cycle.
[nodeListInd1 nodeListInd2] = edges2nodes(edgeListInd_1,:);
nextClockwiseEdgeListInd1 = getNextClockwiseEdge();
% One of the two edges belong to the current region. Keep this edge as the
% next edge

% From this edge, pick the node at the other end find the edge attached to
% it in the same region. This is the next edge. Continue finding the next
% edges in the clockwise order until all the edges are collected.

function nextCWedgeLInd = getNextClockwiseEdges()
