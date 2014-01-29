function ws_grid = getWSfromGrid_sq(nodeIndsVect,edges2pixels,sizeR,sizeC,...
                    gridResolution)

% generates a ws like output from the square grid representation
% all edges and nodes will be labeled zero where as all other 'regions'
% will be given a unique ID, including the image border (=1)

% init
ws_grid = ones(sizeR,sizeC);
% border already assinged label '1'

% mark all edges and nodes with label '0'
ws_grid(nodeIndsVect) = 0;
% remove 1st col (edgeIDs)
edges2pixels(:,1) = [];
% get all edge pixel inds (nonzero)
allEdgePixelInds = edges2pixels(edges2pixels>0);
ws_grid(allEdgePixelInds) = 0;

% assign unique region IDs (starting from 2)
% (using the fact that the grid cells are squares)
