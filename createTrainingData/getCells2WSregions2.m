function [c_wsIDsInCell,c_internalEdgeLIDsInCell,c_extEdgeLIDsInCell,...
          c_internalNodeListInds,c_extNodeListInds]...
            = getCells2WSregions2(labelImg_indexed,ws,numLabels,setOfRegions,...
            edgeListInds,edges2nodes)

% version 2. with directed edges. outputs edgeLIDs instead of eIDs.
        
% Inputs:
%   setOfRegions: contains edgeIDs for each wsRegion
%   nodeEdges: edgeIDs for each nodePixelInd
%   connectedJunctionIDs: 

% Outputs:
%     c_wsIDsInCell
%     c_internalEdgeIDsInCell
%     c_extEdgeIDsInCell
%     c_internalNodeListInds
%     c_extNodeListInds


% Params:
threshFrac = 0.33; % fraction of WSregion pixels allowed to be OUTSIDE of the cell

c_wsIDsInCell = cell(numLabels,1);
c_internalEdgeLIDsInCell = cell(numLabels,1); % off edges - both directions
c_extEdgeLIDsInCell = cell(numLabels,1); % boundary edges LIDs
c_internalNodeListInds = cell(numLabels,1); % off nodes
c_extNodeListInds = cell(numLabels,1);

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
    c_internalEdgeLIDsInCell{i} = internalEdgeIDs_i;
    
    % the other edges are external
    extEdgeIDs_i = setdiff(edgeIDs_unique_i,internalEdgeIDs_i);
    extEdgeLIDs_i = getCwComponentOfEdges(extEdgeIDs_i);
    c_extEdgeLIDsInCell{i} = extEdgeLIDs_i; 
    
    % all nodes where 2 internal edges meet are internal nodes    
    % nodes where at least 1 external edge meets with other nodes, are
    % external nodes
    [internalNodeListInds,extNodeListInds] = getNodes...
                (internalEdgeIDs_i,extEdgeIDs_i,edges2nodes,edgeListInds);
    c_internalNodeListInds{i} = internalNodeListInds;
    c_extNodeListInds{i} = extNodeListInds;
            
end

function extEdgeLIDs_i = getCwComponentOfEdges(extEdgeIDs_i)
% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [internalNodeListInds,extNodeListInds] = getNodes...
                (internalEdgeIDs,extEdgeIDs,edges2nodes,edgeListInds)

intEdgeListInd_logical = ismember(edgeListInds,internalEdgeIDs);            
internalNodeListInds_all = edges2nodes(intEdgeListInd_logical,:);
internalNodeListInds_all = unique(internalNodeListInds_all);

extEdgeListInd_logical = ismember(edgeListInds,extEdgeIDs);            
extNodeListInds = edges2nodes(extEdgeListInd_logical,:);
extNodeListInds = unique(extNodeListInds);

internalNodeListInds = setdiff(internalNodeListInds_all,extNodeListInds);

