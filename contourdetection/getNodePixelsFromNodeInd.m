function nodePix = getNodePixelsFromNodeInd(nodeListInds,nodeIndsAll,connectedJunctionIDsAll)

nodePixInds = nodeIndsAll(nodeListInds);

clusteredNodeInds_logical = ismember(connectedJunctionIDsAll(:,1),nodePixInds);

clusIDs = connectedJunctionIDsAll(clusteredNodeInds_logical,2);

clusNodeInds_logical = ismember(connectedJunctionIDsAll(:,2),clusIDs);

nodePixClustered = connectedJunctionIDsAll(clusNodeInds_logical,1);

nodePix = [nodePixInds; nodePixClustered];

nodePix = unique(nodePix);

