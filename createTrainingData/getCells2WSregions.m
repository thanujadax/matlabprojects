function [c_wsIDsInCell,c_internalEdgeIDsInCell,c_extEdgeIDsInCell]...
            = getCells2WSregions(labelImg_indexed,ws,numLabels,setOfRegions)
        
% Inputs:
%   setOfRegions: contains edgeIDs for each wsRegion

c_wsIDsInCell = cell(numLabels,1);
c_internalEdgeIDsInCell = cell(numLabels,1);
c_extEdgeIDsInCell = cell(numLabels,1);

% regionID = wsID - 1;


parfor i=1:numLabels
    % label_i corresponds to cell_i
    ws_i = ws;
    setOfRegions_i = setOfRegions;
    pixelsForLabel_i = (labelImg_indexed==i);
    wsIDsForLabel_i = ws_i(pixelsForLabel_i);
    wsIDsForLabel_i = unique(wsIDsForLabel_i);
    wsIDsForLabel_i = wsIDsForLabel_i(wsIDsForLabel_i>0);
    c_wsIDsInCell{i} = wsIDsForLabel_i; 
    regionIDs = wsIDsForLabel_i - 1;
    regionIDs = regionIDs(regionIDs>0);
    % get all edges for the cell
    edgesForCell_i = setOfRegions_i(regionIDs,:);
    edgesForCell_i = edgesForCell_i(edgesForCell_i>0);
    
    % pick the internal ones, having 2 regions bounded by it, which are
    % part of this cell
    edgeIDs_unique_i = unique(edgesForCell_i);
    edgeCounts_i = histc(edgesForCell_i,edgeIDs_unique_i);
    internalEdgeIDs_i = edgeIDs_unique_i(edgeCounts_i>1);
    c_internalEdgeIDsInCell{i} = internalEdgeIDs_i;
    
    % the other edges are external
    extEdgeIDs_i = setdiff(edgeIDs_unique_i,internalEdgeIDs_i);
    c_extEdgeIDsInCell{i} = extEdgeIDs_i;
end


