function nextCwEdgeLID = nextCwEdgeTieBreaker(...
    nodeInd,prevEdgeLID,nextEdgeLIDsList,edgepixels,...
    connectedJunctionIDs,sizeR,sizeC)
% when there are multiple edges as the nextCwEdge, it's because
%   1. it's a cluster node with subclusters
%   2. each subcluster will have a legitimate nextCwEdge
% Therefore, we need to pick the edge that is attached to the same
% subcluster as the previous edge

% Inputs:
%   connectedJunctionIDs : first col is pixInds of nodes. second col is
%   clusNodeID for the corresponding pixel
%% identify subclusters for this cluster node
clusID = connectedJunctionIDs((connectedJunctionIDs(:,1)==nodeInd),2);
clusPixInds = connectedJunctionIDs((connectedJunctionIDs(:,2)==clusID),1);
subClusterIDs = getSubClusters(clusPixInds,sizeR,sizeC);
%% get the edges for the subcluster
clusIDsAll = subClusterIDs(:,2);
clusIDsAll = clusIDsAll(clusIDsAll>0);
if(numel(unique(clusIDsAll))==1)
    error('no subclusters found in clustered node for tie breaking nextCwEdgeLIDs!!')
end
% which cluster pixel is the incoming edge connected to
prevEdgePixels = edgepixels(prevEdgeLID,:);
prevEdgePixels = prevEdgePixels(prevEdgePixels>0);

thisClusPixInd = getNodePixForEdge(clusPixInds,prevEdgePixels,sizeR,sizeC);

% get the sub cluster for all the next edges (to be tiebroken)
subClusIDForNextEdges = zeros(numel(nextEdgeLIDsList),1);
for i=1:numel(nextEdgeLIDsList)
    
end


% get the edge connected to the same subcluster
