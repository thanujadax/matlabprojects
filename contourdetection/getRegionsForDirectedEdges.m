function dEdges2regionsOnOff = getRegionsForDirectedEdges...
            (c_edgeLIDsForRegions_cw,edges2nodes_directional,setOfRegions_edgeLIDs)

% corresponding on/off regions for the directional edges in edges2nodes
% the first direction is N1->N2
% dEdges2regionsOnOff = edgeListInd (rowID) | onRegion | offRegion  : dir N1->N2

numRegions = numel(c_edgeLIDsForRegions_cw);
numEdgesDirectional = size(edges2nodes_directional,1);
dEdges2regionsOnOff = zeros(numEdgesDirectional,2);

for i=1:numRegions
    edgeLIDs_region = c_edgeLIDsForRegions_cw{i};
    


end
