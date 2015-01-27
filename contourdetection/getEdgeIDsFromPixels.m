function edgeListInds = getEdgeIDsFromPixels(someEdgePixels,edges2pixels)

% Inputs
%   someEdgePixels - unordered set of edge pixels possibly belonging to
%   multiple edges
%   edges2pixels

% Output
%   edgeLID - edgeListInds

edgepixels = edges2pixels;
edgepixels(:,1) = [];

[~,i_edgepixels,~] = intersect(edgepixels,someEdgePixels);

[sizer,sizec] = size(edgepixels);
[edgeListInds,~] = ind2sub([sizer sizec],i_edgepixels);

edgeListInds = unique(edgeListInds);
