function cellPriors = getCellPriors_intensity(imgIn,setOfCells,edges2pixels,...
    nodeInds,edges2nodes)

[sizeR,sizeC] = size(imgIn);

numCells = size(setOfCells,1);
cellPriors = zeros(numCells,1);

for i=1:numCells
    
    edgeSet_cell = setOfCells(i,:);
    edgeSet_cell = edgeSet_cell(edgeSet_cell>0);
    % get boundary pixels of each cell
    boundaryPixels = getBoundaryPixelsForCell(edgeSet_cell,edges2pixels,...
        nodeInds,edges2nodes,edges2pixels(:,1));

    % get internal pixels of each cell
    [internalx,internaly] = getInternelPixelsFromBoundary(boundaryPixels,sizeR,sizeC);
    
    intPixInds = sub2ind([sizeR sizeC],internaly,internalx);
    
    pixelValues = imgIn(intPixInds);
    
    cellPriors(i) = mean(pixelValues);
    
end

cellPriors = 1 - cellPriors;