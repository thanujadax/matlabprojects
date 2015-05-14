function subClustIDForNextEdge = getSubClustIDForEdge...
    (edgePixInds,clusPixInds,sizeR,sizeC,subClusterIDs)

% Inputs

%   subClusterIDs - temporary integers not to be confused with global
% variable clusterIDs. first column contains the pixelInd. 2nd col contains
% the corresponding temporary subClustID


% Outputs
%   subClustIDForNextEdge 

% To which nodePixInd is the edge attached
nodePixIndForEdge = getNodePixForEdge(clusPixInds,edgePixInds,sizeR,sizeC);

% which subClustID does this nodePixInd have
rowIndLogical = (subClusterIDs(:,1)==nodePixIndForEdge);
subClustIDForNextEdge = 0 ; % init
if(sum(rowIndLogical==1))
    subClustIDForNextEdge = subClusterIDs(rowIndLogical,2);
else
    error('Error looking for matching subclusterID for given edge')
end
