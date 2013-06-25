function setOfCellsMat = setOfCells2Mat(setOfCells)
% takes in a cell array. each cell contains a row vector with some number
% of elements.
% the output is a vertical concatanation of all such row vectors in the
% cell array, padded with zeros at the end

% get the max number of elements in the entire set of row vectors
numCells = numel(setOfCells);
maxNumElements = 0;
for i=1:numCells
    numElementsInCell = numel(setOfCells{i});
    if(maxNumElements<numElementsInCell)
        maxNumElements = numElementsInCell;
    end
end

setOfCellsMat = zeros(numCells,maxNumElements);

% fill in the matrix
for i=1:numCells
    numElementsInCell = numel(setOfCells{i});
    rowVectorFromCell = setOfCells{i};
    for j=1:numElementsInCell
        setOfCellsMat(i,j) = rowVectorFromCell(j);
    end
end