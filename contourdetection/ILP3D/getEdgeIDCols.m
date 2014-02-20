function edgeID_cols_i = getEdgeIDCols(edgeIDs_i,sect1IDs_i,sectVarStats)

% Inputs:
%   sectVarStats - contains the number of edge, node, region vars for each
%       section. rowID=sectID. cols: {edgeVars}{nodeConfs}{regionVars}

offset = 0;

for i=1:(sect1IDs_i-1)
    offset = offset + sum(sectVarStats(i,:));
end

edgeID_cols_i = edgeIDs_i + offset;