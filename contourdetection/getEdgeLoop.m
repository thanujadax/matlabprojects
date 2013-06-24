function setOfEdges = getEdgeLoop(startNodeListInd,startEdgeID,nodeEdges,...
    junctionTypeListInds,jAnglesAll_alpha)
% returns the set of edges that completes a clockwise loop
setOfEdges = [];
nextEdgeID = [];

while(nextEdgeID~=startEdgeID)
    nextEdgeID = [];
    nextEdgeID = getNextEdge(currentEdge,currentNode,nodeEdges,junctionTypeListInds,...
    jAnglesAll_alpha);
    if(~isempty(nextEdgeID))
        setOfEdges = [setOfEdges; nextEdgeID];
    else
        % no_next_edge
        setOfEdges = 0;
        return
    end
end