function [setOfEdges,edgeUsage] = getEdgeLoop(startNodeListInd,startEdgeID,nodeEdges,...
    junctionTypeListInds,jAnglesAll_alpha,edgeUsage,boundaryEdgeIDs,maxUsage,maxUsageB)
% returns the set of edges that completes a clockwise loop
% returns empty if no other edges are found or returns zero if the found
% edges have reached allowed usage limit
% setOfEdges - row vector containing the edgeIDs corresponding to a cell

setOfEdges = [];
nextEdgeID = [];

while(nextEdgeID~=startEdgeID)
    nextEdgeID = getNextEdge(startEdgeID,startNodeListInd,nodeEdges,junctionTypeListInds,...
    jAnglesAll_alpha);
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
    nextEdgeID = [];
end