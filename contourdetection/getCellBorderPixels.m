function cellBorderPixels = getCellBorderPixels(c_cellBorderEdgeIDs,...
            c_cellBorderNodeIDs,edgepixels,nodeInds,connectedJunctionIDs)

numCells = numel(c_cellBorderEdgeIDs);
cellBorderPixels = [];

for i=1:numCells
    % get edge pixels
    cellEdgeInds = c_cellBorderEdgeIDs{i};
    cellBorderEdgePixels_i = edgepixels(cellEdgeInds,:);
    cellBorderEdgePixels_i = cellBorderEdgePixels_i(cellBorderEdgePixels_i>0);
    % cellBorderEdgePixels{i} = cellBorderEdgePixels_i;
    cellBorderPixels = [cellBorderPixels; cellBorderEdgePixels_i];
    
    % get node pixels
    cellNodeInds = c_cellBorderNodeIDs{i};
    cellBorderNodePixels_i = nodeInds(cellNodeInds);
    
    cellBorderPixels = [cellBorderPixels; cellBorderNodePixels_i];
    
    clusIDs = [];
    
    % if this is a cluster node, add the entire cluster to the list
    if(~isempty(connectedJunctionIDs))
        clusID_indList = ismember(connectedJunctionIDs(:,1),cellBorderNodePixels_i);
        clusIDs = connectedJunctionIDs(clusID_indList,2);
        clusIDs = unique(clusIDs);
        
        if(~isempty(clusIDs))
        % cluster node
        clusPixListInds = ismember(connectedJunctionIDs(:,2),clusIDs);
        clusPix_i = connectedJunctionIDs(clusPixListInds,1);
        cellBorderPixels = [cellBorderPixels; clusPix_i];
        end
    end
    
end

cellBorderPixels = unique(cellBorderPixels);