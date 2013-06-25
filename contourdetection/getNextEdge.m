function [nextEdgeId,nextNodeListInd] = getNextEdge(currentEdge,currentNode,nodeEdges,...
    junctionTypeListInds,jAnglesAll_alpha,edges2nodes,edgeIDs)

% Inputs:
%   currentNode - nodeListInd of the current node
%   currentEdge - edgeID of the current edge
%   clockwiseLoop - if 1, we are looking for clockwise loops. this affects
%   the extraction of the next edge at thisNode coming from thisEdge

% To follow the clockwise loop through the current node, we should pick the
% next edge that makes the smallest angle in the counter-clockwise
% direction

% get the edges connected to this node and the angles
% angles
[jListInd,jType] = find(junctionTypeListInds==currentNode);
nodeAnglesAll_alpha = jAnglesAll_alpha{jType};
connectedEdgeAngles_i = nodeAnglesAll_alpha(jListInd,:);
% edge IDs
connectedEdgeIDs_i = nodeEdges(currentNode,:);
connectedEdgeIDs_i = connectedEdgeIDs_i(connectedEdgeIDs_i>0);
connectedEdgeIDs_i(1) = []; % first element is the pixel index of the node

% calculate the angle differences from current edge to other edges
% edgePosThisNode = find(connectedEdgeIDs_i==currentEdge);
% angle_thisEdge = connectedEdgeAngles_i(edgePosThisNode);

angle_thisEdge = connectedEdgeAngles_i(connectedEdgeIDs_i==currentEdge);

angleDiffVector = connectedEdgeAngles_i - angle_thisEdge;
% fixing required since alpha ranges between 0 and 360
negPos = find(angleDiffVector<0);
if(~isempty(negPos))  % if there are negative differences add 360 to those
    numNegPos = numel(negPos);
    for i = 1:numNegPos
        angleDiffVector(negPos(i)) = angleDiffVector(negPos(i)) + 360; 
    end
end
% the edge corresponding to the largest angle difference is the next edge
[~,maxPos] = max(angleDiffVector);
nextEdgeId = connectedEdgeIDs_i(maxPos);
% get the relevant node for the next edge
% need the edgeListInd to get the proper nodeListInd
edgeListInd = find(edgeIDs==nextEdgeId);
nodeInds_nextEdge = edges2nodes(edgeListInd,:);
currentNode_pos = find(nodeInds_nextEdge==currentNode); % is either 1 or 2
if(currentNode_pos==1)
    nextNode_pos = 2;
else
    nextNode_pos = 1;
end
nextNodeListInd = nodeInds_nextEdge(nextNode_pos);