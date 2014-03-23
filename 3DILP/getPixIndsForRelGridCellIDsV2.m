function onCellPixels = getPixIndsForRelGridCellIDsV2...
                        (gridCellRelIDs,rootPixels,gridResY,gridResX,...
                        sizeR,sizeC)

% Input
%   gridCellRelIDs

gridSize = gridResY * gridResX;

numOnCells = numel(gridCellRelIDs);

onCellPixels = zeros(numOnCells,gridSize);

% get root pixels
rootPixForOnCells = rootPixels(gridCellRelIDs);

for i=1:numOnCells
    % get r,c for the entire cell
    [r0,c0] = ind2sub([sizeR sizeC],rootPixForOnCells(i));

    rEnd = r0 + gridResY-1;
    cEnd = c0 + gridResX-1;

    image = zeros(sizeR,sizeC);
    image(r0:rEnd,c0:cEnd) = 1;
    onCellPixels(i,1:gridSize) = find(image);

end