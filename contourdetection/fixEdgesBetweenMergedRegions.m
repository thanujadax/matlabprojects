function [newWS, invisibleEdgeLIDs] = fixEdgesBetweenMergedRegions(...
    newWS,cell_mergedWsIDs_original,wsOld,edges2pixels,expandedWsIDs)

% Inputs:
%   newWS - 
%   mergedWsIDs_original - cell array. each array element is a vector with
%   the old wsRegionIDs
%   wsOld - 
%   edges2pixels - first column contains the edgeIDs
%   expandedWsIDs - the regions that spread to other regions. The order
%   directly corresponds to the merged regions given in
%   cell_mergedWsIDs_original

% Outputs - 

% find out whether there are edges separating merged regions. If so, remove
% them and assign them the mergedRegionID

invisibleEdgeLIDs = [];
edgepixels = edges2pixels;
edgepixels(:,1) = [];

for i = 1:numel(cell_mergedWsIDs_original)
    mergedRegionIDs_orig = cell_mergedWsIDs_original{i};
    
    edgeLIDssCommonToMergedRegions = getCommonEdgesForRegions(...
                mergedRegionIDs_orig,wsOld,edges2pixels);
    
    % fix the edges if any found
    numCommonEdges = sum(edgeLIDssCommonToMergedRegions>0);
    if(numCommonEdges>0)
        invisibleEdgeLIDs = [invisibleEdgeLIDs; edgeLIDssCommonToMergedRegions];
        % fix newWS
        % get edge pixels
        edgePixelsToChange = edgepixels(edgeLIDssCommonToMergedRegions,:);
        edgePixelsToChange = edgePixelsToChange(edgePixelsToChange>0);
        % change ws id
        wsIDToAssign = expandedWsIDs(i);
        newWS(edgePixelsToChange) = wsIDToAssign;
    end
                 
end