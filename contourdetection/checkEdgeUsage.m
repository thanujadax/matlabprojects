function edgeOk = checkEdgeUsage(currentEdgeID,edgeUsage,boundaryEdgeIDs,maxUsage,maxUsageB)
% returns 1 if the edgeUsage is OK. otherwise 0.

edgeListInd = find(edgeUsage(:,1)==currentEdgeID);

if(max(boundaryEdgeIDs==currentEdgeID))
    % currentEdge is a boundary edge
    if(edgeUsage(edgeListInd,2)<maxUsageB)
        edgeOk = 1;
    else
        edgeOk = 0;
    end
else
    % currentEdge is not a boundary edge (most likely)
    if(edgeUsage(edgeListInd,2)<maxUsage)
        edgeOk = 1;
        % flag set. get loop
    else
        edgeOk = 0;
    end
end