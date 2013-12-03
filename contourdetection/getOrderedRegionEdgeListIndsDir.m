function [c_edgeLIDsForRegions_dir_cw,setOfRegions_edgeLIDs] ...
        = getOrderedRegionEdgeListIndsDir...
        (setOfRegions,edges2nodes,jAnglesAll_alpha,...
        junctionTypeListInds,nodeEdgeIDs,edgeListIndsAll)

% Inputs:
%   setOfRegions: edgeIDs for each region as a row vector

% Outputs:
%   c_edgeLIDsForRegions_cw: set of cells each containing directional
%   edgeLIDs for the corresponding region
%   setOfRegions_edgeLIDs: edgeLIDs for each region (undirected)

% at each node take a right turn. get the edge. check if that edge is in
% the set of edges

numRegions = size(setOfRegions,1);
c_edgeLIDsForRegions_dir_cw = cell(numRegions,1);
setOfRegions_edgeLIDs = setOfRegions;

edges2nodes_complements = edges2nodes;
edges2nodes_complements(:,1) = edges2nodes(:,2);
edges2nodes_complements(:,2) = edges2nodes(:,1);
edges2nodes_directional = [edges2nodes; edges2nodes_complements];

for i=1:numRegions
    edgeIDsOfRegion_i = setOfRegions(i,:);
    edgeIDsOfRegion_i = edgeIDsOfRegion_i(edgeIDsOfRegion_i>0);
    % arrange the edges in clockwise order of node traversal
    [~,edgeLIDsForRegion] = intersect(edgeListIndsAll,edgeIDsOfRegion_i);
    
    numE_region = numel(edgeLIDsForRegion);
    setOfRegions_edgeLIDs(i,1:numE_region) = edgeLIDsForRegion;

    nodeLIdsForRegion = edges2nodes_directional(edgeLIDsForRegion,:);
    nodeLIdsForRegion = unique(nodeLIdsForRegion);
    c_edgeLIDsForRegions_dir_cw{i} = getCwOrderedEdgesForRegion...
            (edgeLIDsForRegion,edges2nodes_directional,junctionTypeListInds,...
            nodeEdgeIDs,jAnglesAll_alpha,edgeListIndsAll,nodeLIdsForRegion);
    
end


function cwOrderedEdgeListInds = getCwOrderedEdgesForRegion...
            (edgeListInds_region,edges2nodes_directional,junctionTypeListInds,...
            nodeEdgeIDs,jAnglesAll_alpha,edgeListIndsAll,nodeLIdsForRegion)
% Pick one edge (1st in the list)
edgeLId_1 = edgeListInds_region(1);
edgeID_1 = edgeListInds(edgLId_1);
% Get the nodes at each end of the edge. At each node get the next edge as
% if to complete a clockwise cycle.
[nodeListInd_1,nodeListInd_2] = edges2nodes_directional(edgeLId_1,:);
nextCwEdgeLId_1 = getNextClockwiseEdge(nodeListInd_1,edgeLId_1,edgeID_1,...
            nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListIndsAll);

nextCwEdgeLId_2 = getNextClockwiseEdge(nodeListInd_2,edgeLId_1,edgeID_1,...
            nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListIndsAll);

% One of the two edges belong to the current region. Keep this edge as the
% next edge
nextCwEdgeLId_inRegion = intersect...
                (edgeListInds_region,[nextCwEdgeLId_1,nextCwEdgeLId_2]);

if(isempty(nextCwEdgeLId_inRegion))
    disp('ERROR1: getOrderedRegionEdgeListIndsDir.m nextCwEdge not found!')
else
    % From this edge, pick the node at the other end find the edge attached to
    % it in the same region. This is the next edge. Continue finding the next
    % edges in the clockwise order until all the edges are collected.
    if(nextCwEdgeLId_inRegion==nextCwEdgeLId_1)
        nodeLId_0 = nodeListInd_1;
    else
        nodeLId_0 = nodeListInd_2;
    end
    
    cwOrderedEdgeListInds = getCwSetOfEdges(edgeLId_1,nextCwEdgeLId_inRegion,...
                    nodeLId_0,edges2nodes_directional,edgeListInds_region,...
                    nodeLIdsForRegion);
    
end

function cwOrderedEdgeListInds = getCwSetOfEdges(edgeLId_1,edgeLId_2,...
                        nodeLId_0,edges2nodes_directional,edgeListInds_region,...
                        nodeLIdsForRegion)
numEdges_region = numel(edgeListInds_region);
cwOrderedEdgeListInds = zeros(numEdges_region);

N1 = nodeLId_0;
% N2 is at the other end of the edge edgeLId_2;
nodeLIds_temp = edges2nodes_directional(edgeLId_2,:);
N2 = setDiff(nodeLIds_temp,N1);

%cwOrderedEdgeListInds(1) = edgeLId_2; % is it in the right direction

for i=1:numEdges_region
    % get all edges connected to N2 as the starting node (N1)
    N1 = N2;
    edgesForN1_logical = (edges2nodes_directional(:,1)==N1);
    allN2sForN1 = edges2nodes_directional(edgesForN1_logical,2);
    % which one of those edges belongs to this region
    N2 = intersect(allN2sForN1,nodeLIdsForRegion);
    % which edge is N1->N2: which row has (N1,N2)?
    nextEdgeLID_dir = strmatch([N1,N2],edges2nodes_directional);
    cwOrderedEdgeListInds(i) = nextEdgeLID_dir;    
end


function nextCwEdgeLInd = getNextClockwiseEdge(nodeLId,edgeLId,edgeID,...
            nodeEdgeIDs,junctionTypeListInds,jAnglesAll_alpha,edgeListInds)
% at the node nodeLId, wrt the edge edgeLId, what is the next edge in cw direction
nodeEdgeIDsAll = nodeEdgeIDs(nodeLId,:);
nodeEdgeIDsAll(1) = []; % first element is the nodePixelInd
nodeEdgeIDsAll = nodeEdgeIDsAll(nodeEdgeIDsAll>0);

[junctionListInd,junctionType] = find(junctionTypeListInds==nodeLId);
alphas_junctionType = jAnglesAll_alpha{junctionType};

nodeAlphas_0 = alphas_junctionType(junctionListInd,:);
numOfUniqueAngles = unique(nodeAlphas_0);
if(numel(numOfUniqueAngles)<numel(nodeAlphas_0))
    % duplicate alphas detected. recalculate using only the immediate edge pixel
    % wrt to node pixel
    nodeAlphas_0 = recalculateAlphas(nodeListInd_0,nodeEdges,...
        edges2pixels,sizeR,sizeC);
end
% inputEdgePos = find(nodeEdgesAll==edgeID);
alpha_inputEdge = nodeAlphas_0(nodeEdgeIDsAll==edgeID);
% clockwise search
% get the edge with next largest alpha
% get all angles larger than alpha_0    
largerAlphas = nodeAlphas_0(nodeAlphas_0>alpha_inputEdge);
if(~isempty(largerAlphas))
    nextAlpha = min(largerAlphas);
else
    % get the smallest alpha
    nextAlpha = min(nodeAlphas_0);
end
nextCwEdgeID = nodeEdgeIDsAll(nodeAlphas_0==nextAlpha);
[~,nextCwEdgeLInd] = intersect(edgeListInds,nextCwEdgeID);


function alphas = recalculateAlphas(nodeListInd,nodeEdges,edges2pixels,sizeR,sizeC)
% alphas wrt to the node pixel using only the immediate edge pixels
% inputs:
%   nodeListInd
%   nodeEdges
%   edges2pixels
MAX_NUM_PIXELS = 1;
nodePixInd = nodeEdges(nodeListInd,1);
nodeEdgeIDs = nodeEdges(nodeListInd,:);
nodeEdgeIDs(1) = []; % the first element is node
nodeEdgeIDs = nodeEdgeIDs(nodeEdgeIDs>0);
numEdges = numel(nodeEdgeIDs);
alphas = zeros(1,numEdges);
[y0,x0] = ind2sub([sizeR,sizeC],nodePixInd);

for i=1:numEdges
    edgePixelInds=edges2pixels((edges2pixels(:,1)==nodeEdgeIDs(i)),:);
    edgePixelInds(1) = []; % first element is the edgeID
    edgePixelInds = edgePixelInds(edgePixelInds>0);
    nodeEdgePixel = getNodeEdgePixel(nodePixInd,edgePixelInds,sizeR,sizeC,...
                                    MAX_NUM_PIXELS);
    [y1,x1] = ind2sub([sizeR,sizeC],nodeEdgePixel);
    y = y1 - y0;
    x = x1 - x0;
    
    alpha_i = atan2d(y,x);
    if(alpha_i<0)
        alpha_i = alpha_i + 360;
    end
    
    alphas(i) = alpha_i;

end % for loop
% end of function recalculateAlphas