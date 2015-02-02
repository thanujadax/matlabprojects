function yes = isClusterNode(nodeInd,connectedJunctionIDs)

% inputs:
   % nodeInd - node pixel ind

yes = sum(connectedJunctionIDs(:,1)==nodeInd);