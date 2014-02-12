function jAnglesAll = getNodeAnglesFromEdgeOrientations...
            (junctionTypeListInds,jEdges,edgeOrientationsList,...
            edgeResponses_signed)
        
% Output:
% jAnglesAll: for each node given in the order of junctionTypeListInds,
% stores the orientations of each edge connected to it. equivalent in structure to
% jEdges.

% Inputs:
% edgeOrientationsList: col1: default orientation (N1->N2)