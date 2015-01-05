function edgeListInds = getEdgeIDsFromPixels(someEdgePixels,edges2pixels,...
                    sizeR,sizeC)

% Inputs
%   someEdgePixels - unordered set of edge pixels possibly belonging to
%   multiple edges
%   edges2pixels

% Output
%   edgeLID - edgeListInds

edgepixels = edges2pixels;
edgepixels(:,1) = [];

[~,i_edgepixels,~] = intersect(edgepixels,someEdgePixels);

[edgeListInds,~] = ind2sub([sizeR sizeC],i_edgepixels);

edgeListInds = unique(edgeListInds);
