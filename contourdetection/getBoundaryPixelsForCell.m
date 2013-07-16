function boundaryPixelInds = getBoundaryPixelsForCell(edgeSet_cell,edges2pixels,...
    nodeInds,edges2nodes,edgeIDsAll)

% returns the pixelInds for the boundary pixels for the given cell
boundaryPixelInds = [];

edgeSet_cell = edgeSet_cell(edgeSet_cell>0);

% get node pix inds
numEdges = numel(edgeSet_cell);
edgeListInds = zeros(numEdges,1);

for i=1:numEdges
    edgeListInds(i) = find(edgeIDsAll==edgeSet_cell(i));
    % add node pixels
    clear nodeListInds_i nodePixelInds_i
    nodeListInds_i = edges2nodes(edgeListInds(i),:);
    nodePixelInds_i = nodeInds(nodeListInds_i);
    boundaryPixelInds = [boundaryPixelInds; nodePixelInds_i];
    
    % add edge pixels
    clear edgePixels_i
    edgePixels_i = edges2pixels(edgeListInds(i),:);
    edgePixels_i(1) = []; % first element is the edgeID
    edgePixels_i = edgePixels_i(edgePixels_i>0);
    boundaryPixelInds = [boundaryPixelInds; edgePixels_i'];
    
end

boundaryPixelInds = unique(boundaryPixelInds);