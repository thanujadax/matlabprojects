function [c_edgeLIDsForRegions_dir_cw,setOfRegions_edgeLIDs,edgeLIDs2nodes_directional,...
    dirEdges2regionsOnOff] ...
        = getOrderedRegionEdgeListIndsDir2...
        (setOfRegions,edges2nodes,jAnglesAll_alpha,...
        junctionTypeListInds,nodeEdgeIDs,edgeListIndsAll,...
        edges2pixels,sizeR,sizeC)
    
% Inputs:
% setOfRegions - edgeIDs for each region as a row vector
% edges2nodes - rowID = edgeLID. each row contains the two nodeLIDs for the
% edge
% jAnglesAll_alpha - 
% junctionTypeListInds - 
% nodeEdgeIDs - 
% edgeListIndsAll - 
% edges2pixels - 
% sizeR - 
% sizeC - 


% Outputs:
% c_edgeLIDsForRegions_dir_cw - set of cells each containing directional
%   edgeLIDs for the corresponding region
% setOfRegions_edgeLIDs - edgeLIDs for each region (undirected)
% edgeLIDs2nodes_directional - 
% dirEdges2regionsOnOff - edgeListInd_dir (=rowID) | onRegion | offRegion  : dir N1->N2
%   regionID = 0 is for the image border.