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


