function cellCentroid = getCellCentroid(wsgraph,edgeIDset,edges2pixels,sizeR,sizeC,...
        edges2nodes,nodeInds)

% Inputs:
%   wsgraph - matrix containing non-zero valued pixels for edges and
%   junctions. zeros otherwise
%   edgeIDset - set of edges defining the current cell
%   edges2pixels - first column shows the actual edgeID. each row has the
%   list of pixels for the corresponding edge

% Outputs:
%   cellCentroid = [x; y] coordinates of the approximate centroid of the
%   cell

% generate random pairs of edges (from the set of edges for this cell)
numEdgesInCell = numel(edgeIDset);
edgeCombinationsVector = nchoosek(1:numEdgesInCell,2); % each row has 2 elements
% which are 2 vector indices to get 2 edges

% get the set of points where their perpendicular bisectors meet
for i=1:size(edgeCombinationsVector,1)
    clear edge1_pixels edge2pixels
    edge1_listInd = find(edges2pixels(:,1)==edgeIDset(edgeCombinationsVector(i,1)));
    edge1_pixels = edges2pixels(edge1_listInd,:);
    edge1_pixels(1) = []; % first element is the edgeID. get rid of it.
    edge1_pixels = edge1_pixels(edge1_pixels>0);
    
    edge2_listInd = find(edges2pixels(:,1)==edgeIDset(edgeCombinationsVector(i,2)));
    edge2_pixels = edges2pixels(edge2_listInd,:);
    edge2_pixels(1) = []; % first element is the edgeID. get rid of it.
    edge2_pixels = edge2_pixels(edge2_pixels>0);
    
    [points(i,1),points(i,2)] = getPerpendicularBisectorIntersection...
                (edge1_pixels,edge2_pixels,sizeR,sizeC,edges2nodes,nodeInds);
end
% get the centroid/median of these points -> (x,y)