function cellStates = getCellStatesForActiveEdge(edgeID,cellIDs,edgeOrientation,...
    sizeR,sizeC,edges2pixels,edgeIDs_all,edges2nodes,nodeInds,cells2edges)
% Inputs:
%   edgeID
%   cellIDs = [cell1ID; cell2ID] vector
%   edgeOrientation
%   cells2edges - for each cell, contains the set of edgeIDs as a row
%   vector

% Outputs:
%   cellStates = [cell1; cell2]. state=1-> interior, state=0->membrane

% get which side each cell (just one cell is enough) is wrt the edge
% 

% cell1
cell1_edgeIDs = cells2edges(cellIDs(1));
cell1_Centroid = getCellCentroid(cell1_edgeIDs,edges2pixels,edgeIDs_all,...
    sizeR,sizeC,edges2nodes,nodeInds);

% check which side
cell1_isInterior = checkIfCellIsInterior(); % todo

% check which side is cell 1. only if there's an ambiguity, check for cell2

% cell2
cell2_edgeIDs = cells2edges(cellIDs(2));
cell2_Centroid = getCellCentroid(cell2_edgeIDs,edges2pixels,edgeIDs_all,...
    sizeR,sizeC,edges2nodes,nodeInds);


