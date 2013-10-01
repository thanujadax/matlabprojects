function cellBorderPixels = getCellBorderPixels(c_cellBorderEdgeIDs,...
            c_cellBorderNodeIDs,edgepixels,nodeInds)

numCells = numel(c_cellBorderEdgeIDs);
cellBorderPixels = [];

for i=1:numCells
    % get edge pixels
    cellEdgeInds = c_cellBorderEdgeIDs{i};
    cellBorderEdgePixels_i = edgepixels(cellEdgeInds,:);
    cellBorderEdgePixels_i = cellBorderEdgePixels_i(cellBorderEdgePixels_i>0);
    % cellBorderEdgePixels{i} = cellBorderEdgePixels_i;
    cellBorderPixels = [cellBorderPixels; cellBorderEdgePixels_i'];
    
    % get node pixels
    cellNodeInds = c_cellBorderNodeIDs{i};
    cellBorderNodePixels_i = nodeInds(cellNodeInds);
    cellBorderPixels = [cellBorderPixels; cellBorderNodePixels_i'];
end

