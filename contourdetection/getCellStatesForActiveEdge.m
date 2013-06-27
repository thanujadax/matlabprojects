function cellStates = getCellStatesForActiveEdge(edgeID,cellIDs,edgeOrientation,...
    sizeR,sizeC,edges2nodes,nodeInds)
% Inputs:
%   edgeID
%   cellIDs = [cell1ID; cell2ID] vector
%   edgeOrientation

% Outputs:
%   cellStates = [cell1; cell2]. state=1-> interior, state=0->membrane

% get which side each cell (just one cell is enough) is wrt the edge

