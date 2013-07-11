function [setOfEdges,edgeUsage] = getEdgeLoop(startNodeListInd,startEdgeID,nodeEdges,...
    junctionTypeListInds,jAnglesAll_alpha,edgeUsage,boundaryEdgeIDs,maxUsage,maxUsageB,...
    edges2nodes,edgeIDs)
% returns the set of edges that completes a clockwise loop
% returns empty if no other edges are found or returns zero if the edges that have been 
% found, have reached allowed usage limit
% setOfEdges - row vector containing the edgeIDs corresponding to a cell

setOfEdges = [];
nextEdgeID = 0;
currentEdgeID = startEdgeID;
currentNodeListInd = startNodeListInd;
while(nextEdgeID~=startEdgeID)
    [nextEdgeID,nextNodeListInd] = getNextEdge(currentEdgeID,currentNodeListInd,nodeEdges,junctionTypeListInds,...
    jAnglesAll_alpha,edges2nodes,edgeIDs);
    if(~isempty(nextEdgeID))
        % check usage: if ok append to setOfEdges
        edgeOk = checkEdgeUsage(nextEdgeID,edgeUsage,boundaryEdgeIDs,maxUsage,maxUsageB);
        if(edgeOk)
            setOfEdges = [setOfEdges nextEdgeID];
            % update usage stats
            edgeUsage = updateEdgeUsage(nextEdgeID,edgeUsage);
        else
            setOfEdges = 0;
            return
        end
    else
        % no_next_edge
        setOfEdges = 0;
        return
    end
    currentEdgeID = nextEdgeID;
    currentNodeListInd = nextNodeListInd;
end