function edgeUsage = updateEdgeUsage(currentEdgeID,edgeUsage)
% increments the usage of the edge given by currentEdgeID
% Inputs:
%   edgeUsage: first col - edgeIDs, second col - corresponding usage

edgeListInd = find(edgeUsage(:,1)==currentEdgeID);
currentUsage = edgeUsage(edgeListInd,2);
updatedUsage = currentUsage + 1;
edgeUsage(edgeListInd,2) = updatedUsage;