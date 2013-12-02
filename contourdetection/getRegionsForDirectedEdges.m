function dirEdges2regionsOnOff = getRegionsForDirectedEdges...
            (c_edgeLIDsForRegions_cw,edges2nodes_directional,...
            twoRegionEdgeLIDs)

% Inputs:
%   twoRegionEdges - rowID=eLID; the two columns contain the corresponding
%   2 regions.


% corresponding on/off regions for the directional edges in edges2nodes
% the first direction is N1->N2

% Output: 
% dEdges2regionsOnOff = edgeListInd_dir (=rowID) | onRegion | offRegion  : dir N1->N2

numRegions = numel(c_edgeLIDsForRegions_cw);
numEdgesDirectional = size(edges2nodes_directional,1);
dirEdges2regionsOnOff = zeros(numEdgesDirectional,2);

% rID=1 is the image border

for i=2:numRegions
    edgeLIDs_dir_region = c_edgeLIDsForRegions_cw{i};
    
    % get logical indices for edgeLID_dir
    numDirEdge_region = numel(edgeLIDs_dir_region);
    edgeLID_dir_sequence = 1:numDirEdge_region;
    edgeLIDs_dir_region_logical = (edgeLID_dir_sequence==edgeLIDs_dir_region);
    
    dirEdges2regionsOnOff(edgeLIDs_dir_region_logical,1) = i; % =rID_on
    
    % get rID_off
    rID_off = getOffRID(twoRegionEdgeLIDs,i);
    dirEdges2regionsOnOff(edgeLIDs_dir_region_logical,2) = rID_off; % rID_off  

end


function rID_off = getOffRID(twoRegionEdgeLIDs,rID_on)
% col1
col1RID_on_logical = (twoRegionEdgeLIDs(:,1)==rID_on);
rID_off1 = twoRegion(col1RID_on_logical,2);

col2RID_on_logical = (twoRegionEdgeLIDs(:,2)==rID_on);
rID_off2 = twoRegion(col2RID_on_logical,1);

rID_off = [rID_off1; rID_off2];


