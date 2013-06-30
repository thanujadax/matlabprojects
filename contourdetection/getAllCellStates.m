function cellStates = getAllCellStates(cells2edges,cellcogs,edgeIDs,edgeOrientations,...
    sizeR,sizeC,edges2pixels)

% Inputs:
%   cells2edges: each row contains the set of edges bounding the cell given
%   by the row number
%   cellcogs: centre of gravities. structure array
%   edgeIDs
%   edgeOrientations
%   sizeR
%   sizeC
%   edges2pixels

% Output:
%   returns a vector containing the cell state (cell interior = 1; membrane
%   = 0) for each cell


numCells = size(cells2edges,1);
cellStates = zeros(numCells,1);

for i=1:numCells
    clear getEdgeSet_cell 
    getEdgeSet_cell = cells2edges(i,:);
    getEdgeSet_cell = getEdgeSet_cell(getEdgeSet_cell>0);
    numEdges_cell = numel(getEdgeSet);
    clear cellStateVector_i
    cellStateVector_i = zeros(numEdges_cell,1);
    for j=1:numEdges_cell
        edgeListInd_j = find(edgeIDs==getEdgeSet_cell(j));
        clear edgePixels
        edgePixels = edges2pixels(edgeListInd_j,:);
        edgePixels = edgePixels(edgePixels>0);         
        cellStateVector_i(j)=checkIfCellIsInterior...
            (cellcogs(i),edgePixels,edgeOrientations(edgeListInd_j),...
            sizeR,sizeC);
    end
    cellState_i = mean(cellStateVector_i);
    if(cellState_i>0)
        cellStates(i) = 1;
    end
end
