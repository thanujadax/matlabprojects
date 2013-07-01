function cellCentroid = getCellCentroid(edgeIDset_cell,edges2pixels,edgeIDs_all,...
    sizeR,sizeC,edges2nodes,nodeInds)

% Inputs:
%   wsgraph - matrix containing non-zero valued pixels for edges and
%   junctions. zeros otherwise
%   edgeIDset - set of edges defining the current cell
%   edges2pixels - first column shows the actual edgeID. each row has the
%   list of pixels for the corresponding edge

% Outputs:
%   cellCentroid = [x; y] coordinates of the approximate centroid of the
%   cell

% get all the edges bounding the cell
numEdgesForCell = numel(edgeIDset_cell);
edgeListInds = zeros(numEdgesForCell,1);
X = [];
Y = [];
for i=1:numEdgesForCell
    edgeListInds(i) = find(edgeIDs_all==edgeIDset_cell(i));
    clear edge_i_pixels
    edge_i_pixels = edges2pixels(edgeListInds(i),:);
    edge_i_pixels(1) = []; % first element is the edgeID
    edge_i_pixels = edge_i_pixels(edge_i_pixels>0);
    clear x y
    [y x] = ind2sub([sizeR sizeC],edge_i_pixels');
    X = [X; x];
    Y = [Y; y];
end

meanx = floor(mean(X));
meany = floor(mean(Y));

cellCentroid = [meanx meany];