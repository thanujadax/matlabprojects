function [edgeID_cols_i,nodeID_cols_i] = getEdgeNodeIDCols...
                (edgeIDs_i,nodeIDs_i,sect1IDs_i,sectVarStats)

% Inputs:
%   sectVarStats - contains the number of edge, node, region vars for each
%       section. rowID=sectID. cols: {edgeVars}{nodeConfs}{regionVars}

offsetForEdgeIDs = 0;

for i=1:(sect1IDs_i-1)
    offsetForEdgeIDs = offsetForEdgeIDs + sum(sectVarStats(i,:));
end

edgeID_cols_i = edgeIDs_i + offsetForEdgeIDs;

nodeID_cols_i = nodeIDs_i + offsetForEdgeIDs + sectVarStats(sect1IDs_i,1);
