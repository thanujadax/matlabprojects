function onCellPixels = getPixIndsForRelGridCellIDs...
                        (gridCellStates,rootPixels,gridResY,gridResX,...
                        sizeR,sizeC)

% Input
%   gridCellStates - the state vector for gridCells. 1 for exterior. 0
%   for interior

gridSize = gridResY * gridResX;

cellInterior_logicalInd = gridCellStates==0;
numOnCells = sum(cellInterior_logicalInd);
numCellInterior = sum(cellInterior_logicalInd);

onCellPixels = zeros(numCellInterior,gridSize);

% get root pixels
rootPixForOnCells = rootPixels(logical(cellInterior_logicalInd));

for i=1:numOnCells
    % get r,c for the entire cell
    [r0,c0] = ind2sub([sizeR sizeC],rootPixForOnCells(i));

    rEnd = r0 + gridResY-1;
    cEnd = c0 + gridResX-1;

    image = zeros(sizeR,sizeC);
    image(r0:rEnd,c0:cEnd) = 1;
    onCellPixels(i,1:gridSize) = find(image);

end
