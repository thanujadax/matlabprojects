function [newWS, removedWsIDs] = getCorrectedWSregions(ws,removeEdgeIDs,edges2pixels)
% merges ws regions after removing removeEdgeIDs

% Inputs
%   ws - ws transform of the original image
%   removeEdgeIDs - edges to be removed from ws to merge regions
%   edges2pixels - first col has edgeID. other cols have pixelInds

% Outputs
%   newWS - new ws transform after merging ws region due to removal of
%   edges
%   removedWsIDs - list of removed ws region ids due to merging with larger regions

newWS = ws;
edgeIDsAll = edges2pixels(:,1);
edgepixelsAll = edges2pixels;
edgepixelsAll(:,1) = [];
[sizeR,sizeC] = size(ws);
% for each removed edge
for i=1:numel(removeEdgeIDs)
% get the two ws regions on either side
    edgeListInd_i = find(edgeIDsAll==removeEdgeIDs(i));
    edgepixels_i = edgepixelsAll(edgeListInd_i,:);
    edgepixels_i = edgepixels_i(edgepixels_i>0);
    regionIDs = getRegionsForEdgePixels(ws,edgepixels_i,sizeR,sizeC);
% find the bigger region and assign its wsID to the smaller region except
% if its region 1 = border (takes precedence)

end