function nextCwEdgeLID = nextCwEdgeTieBreaker(...
    nodeInd,prevEdgeLID,nextEdgeLIDsList,edgepixels,...
    connectedJunctionIDs,sizeR,sizeC)
% when there are multiple edges as the nextCwEdge, it's because
%   1. it's a cluster node with subclusters
%   2. each subcluster will have a legitimate nextCwEdge
% Therefore, we need to pick the edge that is attached to the same
% subcluster as the previous edge (? are we sure about this ???)

% Inputs:
%   connectedJunctionIDs : first col is pixInds of nodes. second col is
%   clusNodeID for the corresponding pixel
%% identify subclusters for this cluster node
clusID = connectedJunctionIDs((connectedJunctionIDs(:,1)==nodeInd),2);
clusPixInds = connectedJunctionIDs((connectedJunctionIDs(:,2)==clusID),1);
% subClusterIDs are temporary integers not to be confused with global
% variable clusterIDs. first column contains the pixelInd. 2nd col contains
% the corresponding temporary subClustID.
subClusterIDs = getSubClusters(clusPixInds,sizeR,sizeC);
%% get the edges for the subcluster
subClusIDsAll = subClusterIDs(:,2);
subClusIDsAll = subClusIDsAll(subClusIDsAll>0);
if(numel(unique(subClusIDsAll))==1)
    error('no subclusters found in clustered node for tie breaking nextCwEdgeLIDs!!')
end
% which cluster pixel is the incoming edge connected to
prevEdgePixels = edgepixels(prevEdgeLID,:);
prevEdgePixels = prevEdgePixels(prevEdgePixels>0);

% nodePixIndToIncomingEdge = getNodePixForEdge(clusPixInds,prevEdgePixels,sizeR,sizeC);

subClustIDForPrevEdge = getSubClustIDForEdge...
            (prevEdgePixels,clusPixInds,sizeR,sizeC,subClusterIDs);

% get the sub cluster for all the next edges (to be tiebroken)
subClustIDForNextEdges = zeros(numel(nextEdgeLIDsList),1);
for i=1:numel(nextEdgeLIDsList)
    edgePixInds_i = edgepixels(nextEdgeLIDsList(i),:);
    edgePixInds_i = edgePixInds_i(edgePixInds_i>0);
    subClustIDForNextEdges(i) = getSubClustIDForEdge...
            (edgePixInds_i,clusPixInds,sizeR,sizeC,subClusterIDs);
    
end

% get the edge connected to the same subcluster
nextCwEdgeLID = nextEdgeLIDsList(subClustIDForNextEdges==subClustIDForPrevEdge);
if(sum(subClustIDForNextEdges==subClustIDForPrevEdge) ~=1)
    error('tie breaking did not work :-( ')
end
