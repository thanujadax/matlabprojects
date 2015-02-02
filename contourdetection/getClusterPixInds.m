function clusterPixInds = getClusterPixInds(nodeInd,connectedJunctionIDs)

% inputs
% nodeInd - pix ind of reference node

% get clus node ind
clusInd = connectedJunctionIDs((connectedJunctionIDs(:,1)==nodeInd),2);

% get all clus pixels
clusterPixInds = connectedJunctionIDs((connectedJunctionIDs(:,2)==clusInd),1);
