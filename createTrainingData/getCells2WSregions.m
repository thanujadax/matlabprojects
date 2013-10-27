function [c_wsIDsInCell,c_internalEdgeIDsInCell,c_extEdgeIDsInCell]...
            = getCells2WSregions(labelImg_indexed,ws,numLabels,setOfRegions)
        
% Inputs:
%   setOfRegions: contains edgeIDs for each wsRegion


% Params:
threshFrac = 0.2; % fraction of WSregion pixels allowed to be OUTSIDE of the cell

c_wsIDsInCell = cell(numLabels,1);
c_internalEdgeIDsInCell = cell(numLabels,1);
c_extEdgeIDsInCell = cell(numLabels,1);

% regionID = wsID - 1;


parfor i=1:numLabels
    % label_i corresponds to cell_i
    ws_i = ws;
    setOfRegions_i = setOfRegions;
    pixelsForLabel_i_logical = (labelImg_indexed==i);
    wsIDsForLabel_i = ws_i(pixelsForLabel_i_logical);
    wsIDsForLabel_i = unique(wsIDsForLabel_i);
    wsIDsForLabel_i = wsIDsForLabel_i(wsIDsForLabel_i>0);
    
    % pick the wsRegions having 4/5 of their pixels in this cell
    wsIDsForLabel_i_filtered = getGoodWSregions...
        (wsIDsForLabel_i,ws,pixelsForLabel_i_logical,threshFrac);
    
    c_wsIDsInCell{i} = wsIDsForLabel_i_filtered; 
    regionIDs = wsIDsForLabel_i_filtered - 1;
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


function wsIDs_filtered = getGoodWSregions(wsIDsForLabel_i,ws,pixForLabel_j_logical,threshFrac)
wsIDs_filtered = [];
numWsRegions = numel(wsIDsForLabel_i);
for i=1:numWsRegions
    wsID_i = wsIDsForLabel_i(i);
    wsPixels_i_logical = (ws==wsID_i);
    wsPixOutsideLabelj_logical = ((pixForLabel_j_logical - wsPixels_i_logical)==-1);
    numWsPixOutside = sum(sum(wsPixOutsideLabelj_logical));
    numWsPixTot = sum(sum(wsPixels_i_logical));
    fractionOutside = numWsPixOutside/numWsPixTot;
    if (fractionOutside<threshFrac)
        wsIDs_filtered = [wsIDs_filtered; wsID_i];
    end
end