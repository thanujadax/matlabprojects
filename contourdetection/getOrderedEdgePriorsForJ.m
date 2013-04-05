function edgePriors_i = getOrderedEdgePriorsForJ(jTypeId,junctionTypeListInds,...
                    nodeEdges,edgePriors,edgeListInds)
% for the junction type given (i) for all the nodes, return an array of the
% edgePriors in the same order for the edges given in nodeEdges

numEdges = jTypeId + 1;
nodeListInds = junctionTypeListInds(:,jTypeId);
nodeListInds = nodeListInds(nodeListInds>0);
numNodes = numel(nodeListInds);

edgePriors_i = zeros(numNodes,numEdges);

for i=1:numNodes
    nodeListInd = nodeListInds(i);      % list index for this node
    nodeEdgeSet = nodeEdges(nodeListInd,2:(numEdges+1)); % edgeIDs for its edges
    % edgeInds = zeros(numEdges,1);       % to store their edgeListInds
    for j=1:numEdges
        edgeListInd_j = find(edgeListInds==nodeEdgeSet(j));
        edgePriors_i(i,j) = edgePriors(edgeListInd_j);
    end
    
end