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

MAX_NUM_EDGES_TO_CONSIDER = 4;

% generate random pairs of edges (from the set of edges for this cell)
numEdgesInCell = numel(edgeIDset_cell);

numEdgesInCell = min(numEdgesInCell,MAX_NUM_EDGES_TO_CONSIDER);
edgeCombinationsVector = nchoosek(1:numEdgesInCell,2); % each row has 2 elements
% which are 2 vector indices to get 2 edges

% get the set of points where their perpendicular bisectors meet
numCombinations = size(edgeCombinationsVector,1);
points = zeros(numCombinations,2);
for i=1:numCombinations
    clear edge1_pixels edge2pixels
    edge1_ind = edgeIDset_cell(edgeCombinationsVector(i,1));
    edge1_listInd = find(edges2pixels(:,1)==edge1_ind);
    edge1_pixels = edges2pixels(edge1_listInd,:);
    edge1_pixels(1) = []; % first element is the edgeID. get rid of it.
    edge1_pixels = edge1_pixels(edge1_pixels>0);
    
    edge2_ind = edgeIDset_cell(edgeCombinationsVector(i,2));
    edge2_listInd = find(edges2pixels(:,1)==edge2_ind);
    edge2_pixels = edges2pixels(edge2_listInd,:);
    edge2_pixels(1) = []; % first element is the edgeID. get rid of it.
    edge2_pixels = edge2_pixels(edge2_pixels>0);
    
    [x,y] = getPerpendicularBisectorIntersection(edge1_pixels,edge2_pixels,...
    edge1_ind,edge2_ind,sizeR,sizeC,edges2nodes,nodeInds,edgeIDs_all);

    points(i,1) = -1;
    points(i,2) = -1;
    if(~isempty(x) && ~isempty(y))
        if(isnumeric(x) && isnumeric(y))
            points(i,1) = x;
            points(i,2) = y;
        end
    end
        
end
% remove any negative points from the points(i,j) matrix
positiveRowInds = points(:,1)>0;
points_pos = points(positiveRowInds,:);
% get the centroid/median of these points -> (x,y)
cellCentroid = floor(median(points_pos));