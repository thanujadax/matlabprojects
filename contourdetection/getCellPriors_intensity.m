function cellPriors = getCellPriors_intensity(imgIn,setOfCells,edges2pixels,...
    nodeInds,edges2nodes,K)
% Inputs:
%   imgIn - normalized image. 1 -> bright
%   K - positive scalar factor for the costs 

% Output:
%   cellPriors - value between K (membrane) and -K (cell interior). 

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
    if(~isempty(intPixInds))
        pixelValues = imgIn(intPixInds);

        % cellPriors(i) = 1 - 2* mean(pixelValues);

        meanPixVal = mean(pixelValues); % always between 0 and 1
        theta = meanPixVal * pi;
        cellPriors(i) = cos(theta) * K; 
    else
        cellPriors(i) = 0;
    end
end


