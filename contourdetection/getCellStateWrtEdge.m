function cellState = getCellStateWrtEdge(edgeID,edgeSet_cell,...
                    edgeOrientation,edges2pixels,...
                    alphas_all,nodeEdges,junctionTypeListInds,edges2nodes,...
                    sizeR,sizeC)
% returns +1 if it's interior. else -1.
                

% get nodes (we need just 1)
edgeListIndsAll = edges2pixels(:,1);
edgeListInd_0 = find(edgeListIndsAll==edgeID);
nodeListIndsForEdge_0 = edges2nodes(edgeListInd_0,:);

% get alphas for the edges connected to the node
nodeListInd_0 = nodeListIndsForEdge_0(1);  % we need to analyze only one node
nodeEdgeIDs_0 = nodeEdges(nodeListInd_0,:);
nodeEdgeIDs_0(1) = []; % first element is node pixel Ind
nodeEdgeIDs_0 = nodeEdgeIDs_0(nodeEdgeIDs_0>0);
[junctionListInd,junctionType] = find(junctionTypeListInds==nodeListInd_0);
alphas_junctionType = alphas_all{junctionType};
nodeAlphas_0 = alphas_junctionType(junctionListInd,:);
numOfUniqueAngles = unique(nodeAlphas_0);
if(numel(numOfUniqueAngles)<numel(nodeAlphas_0))
    % duplicate alphas detected. recalculate using only the immediate edge pixel
    % wrt to node pixel
    nodeAlphas_0 = recalculateAlphas(nodeListInd_0,nodeEdges,...
        edges2pixels,sizeR,sizeC);
end

% get direction to search (clockwise or counter clockwise)
% to get the next edge which is part of an interior.
% if that is part of the current cell, then the current cell is interior
inputEdgePos = find(nodeEdgeIDs_0==edgeID);
alpha_inputEdge = nodeAlphas_0(inputEdgePos);
directionToSearch = cosd(edgeOrientation - alpha_inputEdge);

% get the cell to the brighter side
% hence the cellState

% get next edge ID
if(directionToSearch>0)
    % clockwise search
    % get the edge with next largest alpha
    % get all angles larger than alpha_0    
    largerAlphas = nodeAlphas_0(nodeAlphas_0>alpha_inputEdge);
    if(~isempty(largerAlphas))
        nextAlpha = min(largerAlphas);
    else
        % get the smallest alpha
        nextAlpha = min(nodeAlphas_0);
    end
else
    % counter-clockwise search
    smallerAlphas = nodeAlphas_0(nodeAlphas_0<alpha_inputEdge);
    if(~isempty(smallerAlphas))
        nextAlpha = max(smallerAlphas);
    else
        % get the largest alpha
        nextAlpha = max(nodeAlphas_0);
    end
end

% which region is cell interior. look for the next edge
nextEdgePos = find(nodeAlphas_0==nextAlpha);
nextEdgeID = nodeEdgeIDs_0(nextEdgePos);
edgeInCell = find(edgeSet_cell==nextEdgeID);
if(~isempty(edgeInCell))
    cellState = 1;
else
    cellState = -1;
end

end % main function

function alphas = recalculateAlphas(nodeListInd,nodeEdges,edges2pixels,sizeR,sizeC)
% alphas wrt to the node pixel using only the immediate edge pixels
% inputs:
%   nodeListInd
%   nodeEdges
%   edges2pixels
MAX_NUM_PIXELS = 1;
nodePixInd = nodeEdges(nodeListInd,1);
nodeEdgeIDs = nodeEdges(nodeListInd,:);
nodeEdgeIDs(1) = []; % the first element is node
nodeEdgeIDs = nodeEdgeIDs(nodeEdgeIDs>0);
numEdges = numel(nodeEdgeIDs);
alphas = zeros(1,numEdges);
[y0,x0] = ind2sub([sizeR,sizeC],nodePixInd);

for i=1:numEdges
    edgePixelInds=edges2pixels((edges2pixels(:,1)==nodeEdgeIDs(i)),:);
    edgePixelInds(1) = []; % first element is the edgeID
    edgePixelInds = edgePixelInds(edgePixelInds>0);
    nodeEdgePixel = getNodeEdgePixel(nodePixInd,edgePixelInds,sizeR,sizeC,...
                                    MAX_NUM_PIXELS);
    [y1,x1] = ind2sub([sizeR,sizeC],nodeEdgePixel);
    y = y1 - y0;
    x = x1 - x0;
    
    alpha_i = atan2d(y,x);
    if(alpha_i<0)
        alpha_i = alpha_i + 360;
    end
    
    alphas(i) = alpha_i;

end % for loop

end % function recalculateAlphas