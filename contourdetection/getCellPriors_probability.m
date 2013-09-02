function cellPriors = getCellPriors_probability(pixelProbabilities,setOfCells,edges2pixels,...
    nodeInds,edges2nodes,K,sizeR,sizeC,wsIndsForCells,ws)
% Inputs:
%   imgIn(pixelProbabilities) - normalized image. 1 -> bright
%   K - positive scalar factor for the costs 

% Output:
%   cellPriors - value between K (membrane) and -K (cell interior). 

% numpixels = numel(pixelProbabilities);

numCells = size(setOfCells,1);
cellPriors = zeros(numCells,1);
regionScoreSpace = zeros(sizeR,sizeC);

for i=1:numCells
    
    edgeSet_cell = setOfCells(i,:);
    edgeSet_cell = edgeSet_cell(edgeSet_cell>0);
    % get boundary pixels of each cell
%     boundaryPixels = getBoundaryPixelsForCell(edgeSet_cell,edges2pixels,...
%         nodeInds,edges2nodes,edges2pixels(:,1));

    % get internal pixels of each cell
%     [internalx,internaly] = getInternelPixelsFromBoundary(boundaryPixels,sizeR,sizeC);
%     
%     intPixInds = sub2ind([sizeR sizeC],internaly,internalx);
    
    intPixInds = getInternalPixForCell(ws,wsIndsForCells(i));
    
    
    if(~isempty(intPixInds))
        pixelValues = pixelProbabilities(intPixInds);

        % cellPriors(i) = 1 - 2* mean(pixelValues);

        meanPixVal = mean(pixelValues); % always between 0 and 1
        theta = meanPixVal * pi;
        cellPriors(i) = cos(theta) * K; 
    else
        cellPriors(i) = 0;
    end
    regionScoreSpace(intPixInds) = cellPriors(i);
end

% visualize region scores
regionScoreSpace = regionScoreSpace./(max(max(regionScoreSpace)));
figure;imshow(regionScoreSpace);title('region scores')
