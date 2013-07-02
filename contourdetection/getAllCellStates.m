function cellStates = getAllCellStates(cells2edges,edgeIDs,edgeOrientations,...
    sizeR,sizeC,edges2pixels,nodeInds,edges2nodes)

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

MAX_NUM_EDGES_TO_CONSIDER = 7;
NUM_INTERNAL_POINTS = 7;

numCells = size(cells2edges,1);
cellStates = zeros(numCells,1);

for i=1:numCells
    clear edgeSet_cell 
    edgeSet_cell = cells2edges(i,:);
    edgeSet_cell = edgeSet_cell(edgeSet_cell>0);
    numEdges_cell = numel(edgeSet_cell);
    clear cellStateVector_i
    cellStateVector_i = zeros(numEdges_cell,1);
    
    numEdges_toConsider = min(numEdges_cell,MAX_NUM_EDGES_TO_CONSIDER);
    % get cellInteriorPoints (10points)
    boundaryPixels = getBoundaryPixelsForCell(edgeSet_cell,edges2pixels,...
    nodeInds,edges2nodes,edgeIDs);
    [internalx,internaly] = getInternelPixelsFromBoundary(boundaryPixels,sizeR,sizeC);
    cellInteriorPoints = samplePointsFromArrays(internalx,internaly,NUM_INTERNAL_POINTS);
    
    for j=1:numEdges_toConsider
        edgeListInd_j = find(edgeIDs==edgeSet_cell(j));
        clear edgePixels
        edgePixels = edges2pixels(edgeListInd_j,:);
        edgePixels = edgePixels(edgePixels>0); 
               
        cellStateVector_i(j)=checkIfCellIsInterior...
            (cellInteriorPoints,edgePixels,edgeOrientations(edgeListInd_j),...
            sizeR,sizeC);
        
    end
    
    cellState_i = mean(cellStateVector_i);
    
    if(cellState_i>0)
        cellStates(i) = 1;
    end
end

end

function cellInteriorPoints = samplePointsFromArrays(internalx,internaly,maxn)
    numInputPoints = numel(internalx);
    % pick maxn number of random integers ranging from 1 to numInputPoints
    pixListInds = randi(numInputPoints,maxn,1);
    cellInteriorPoints = zeros(maxn,2);
    for i=1:maxn
        cellInteriorPoints(i,1) = internalx(pixListInds(i));
        cellInteriorPoints(i,2) = internaly(pixListInds(i));
    end
end

