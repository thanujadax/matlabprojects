function cellState = getCellStateWrtEdge(edgeID,edgeSet_cell,...
                    edgeOrientation,sizeR,sizeC,edges2pixels,...
                    nodeInds,edges2nodes,alphas_all,nodeEdges,junctionTypeListInds)
                

% get nodes (we need just 1)
edgeListIndsAll = edges2pixels(:,1);
edgeListInd_0 = find(edgeListIndsAll==edgeID);
nodeIndForEdge_0 = nodeInds(edgeListInd_0,:);

% get alphas for the edges connected to the node
nodeListInd_0 = nodeInds(1);  % we need to analyze only one node
nodeEdgeIDs_0 = nodeEdges(nodeEdgeIDs_0,:);
nodeEdgeIDs_0(1) = [];
nodeEdgeIDs_0 = nodeEdgeIDs_0(nodeEdgeIDs_0>0);
[junctionListInd,junctionType] = find(junctionTypeListInds==nodeListInd_0);
alphas_junctionType = alphas_all{junctionType};
nodeAlphas_0 = alphas_junctionType(junctionListInd,:);
nodeAlphas_0 = nodeAlphas_0(nodeAlphas_0>0);
% get dirction to search (clockwise or counter clockwise)
inputEdgePos = find(nodeEdgeIDs_0==edgeID);
alpha_inputEdge = nodeAlphas_0(inputEdgePos);
directionToSearch = cosd(edgeOrientation - alpha_inputEdge);

% get the cell to the brighter side
% hence the cellState
numEdgesForNode = numel(nodeEdgeIDs_0);

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
    cellState = 0;
end