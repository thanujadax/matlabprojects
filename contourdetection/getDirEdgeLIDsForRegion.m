function dirEdgeLIDsForRegion = getDirEdgeLIDsForRegion(...
            regionID,startEdgeLID,setOfRegions,edgeLIDs2nodes_directional,...
            edgeListIndsAll)
        
% Inputs:
% regionID - 
% startEdgeLID - 
% setOfRegions - edgeIDs for each region

% Output:

if(regionID==0)
    error('method not defined for the 0 region - image border')
end

% get all edgeLIDs for this region
edgeIDsRegion = setOfRegions(regionID,:);
edgeIDsRegion = edgeIDsRegion(edgeIDsRegion>0);
[~,edgeLIDsRegion] = intersect(edgeListIndsAll,edgeIDsRegion);
numEdgesRegion = numel(edgeLIDsRegion);
dirEdgeLIDsForRegion = zeros(1,numEdgesRegion);

% get the nodes for all edgeLIDs
regionNodeLIDs = edgeLIDs2nodes_directional(edgeLIDsRegion,:);

% get the nodes for the startEdgeLID
nextEdgeNodes = edgeLIDs2nodes_directional(startEdgeLID,:);

% traverse the list of nodeLIDs in the direction of the startEdgeLID
% get the corresponding directedEdgeLIDs

nextEdgeLID = startEdgeLID;

for i=1:numEdgesRegion
    
    % the 2nd node of the prevEdge is the 1st node of the nextEdge
    node1 = nextEdgeNodes(2);
    [~,nextEdgeLIDs] = intersect(edgeLIDs2nodes_directional(:,1),node1);
    dirEdgeLIDsForRegion(i) = nextEdgeLID;
    
end

