function nextEdgeId = getNextEdge(currentEdge,currentNode,nodeEdges,junctionTypeListInds,...
    jAnglesAll_alpha,clockwiseLoop)

% Inputs:
%   currentNode - nodeListInd of the current node
%   clockwiseLoop - if 1, we are looking for clockwise loops. this affects
%   the extraction of the next edge at thisNode coming from thisEdge

% To follow the clockwise loop through the current node, we should pick the
% next edge that makes the smallest angle in the counter-clockwise
% direction

% get the edges connected to this node and the angles
[~,jType] = find(junctionTypeListInds,currentNode);
nodeAngles_alpha = junctionTypeListInds(jType);
connectedEdgeIDs = nodeEdges(currentNode,:);
connectedEdgeIDs = connectedEdgeIDs(connectedEdgeIDs>0);
connectedEdgeIDs(1) = []; % first element is the pixel index of the node

% calculate the angle differences from current edge to other edges
