function jAnglesAll = getNodeAnglesFromEdgeOrientations...
            (jEdges,edgeOrientationsList)
        
% Output:
% jAnglesAll: for each node given in the order of junctionTypeListInds,
% stores the orientations of each edge connected to it. equivalent in structure to
% jEdges.

% Inputs:
% edgeOrientationsList: orientations for positive edge response
% edgeResponses: 
% jEdges: edgeIDs at each node. Organized in a cell array. Each cell
% corresponds to a jType (numEdgesPerJn-1)

[~,numJtypes]=size(jEdges);
jAnglesAll = {};

for i=1:numJtypes
    jEdges_i = jEdges{i};
    
    jAngles_i = edgeOrientationsList(jEdges_i);
    jAnglesAll{i} = jAngles_i;
end
