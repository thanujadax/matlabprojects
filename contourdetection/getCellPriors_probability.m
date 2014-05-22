function regionPriors = getCellPriors_probability(pixelProbabilities,setOfCells,...
    sizeR,sizeC,wsIndsForRegion,ws,displayImg)
% Inputs:
%   imgIn(pixelProbabilities) - normalized image. 1 -> bright
%   K - positive scalar factor for the costs 

% Output:
%   cellPriors - value between K (membrane) and -K (cell interior). 

% numpixels = numel(pixelProbabilities);

numCells = size(setOfCells,1);
regionPriors = zeros(numCells,1);
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
    
    intPixInds = getInternalPixForCell(ws,wsIndsForRegion(i));
    
    
    if(~isempty(intPixInds))
        pixelValues = pixelProbabilities(intPixInds);

        % cellPriors(i) = 1 - 2* mean(pixelValues);

        meanPixVal = mean(pixelValues); % always between 0 and 1
        theta = meanPixVal * pi;
        % regionPriors(i) = cos(theta) * K;
        regionPriors(i) = meanPixVal;
    else
        regionPriors(i) = 0;
    end
    regionScoreSpace(intPixInds) = regionPriors(i);
end

% visualize region scores
regionScoreSpace = regionScoreSpace./(max(max(regionScoreSpace)));
if(displayImg)
    figure;imshow(regionScoreSpace);title('region scores')
end