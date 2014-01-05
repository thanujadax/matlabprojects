function dirEdges2regionsOnOff = getRegionsForDirectedEdges...
            (c_edgeLIDsForRegions_cw,edges2nodes_directional,...
            setOfRegions_edgeLIDs,numEdges)

% Inputs:
%   twoRegionEdges - rowID=eLID; the two columns contain the corresponding
%   2 regions.
%   edgeIDs2regions


% corresponding on/off regions for the directional edges in edges2nodes
% the first direction is N1->N2

% Output: 
% dEdges2regionsOnOff = edgeListInd_dir (=rowID) | onRegion | offRegion  : dir N1->N2
%   regionID = 0 is for the image border.

numRegions = numel(c_edgeLIDsForRegions_cw);
numEdgesDirectional = size(edges2nodes_directional,1);
dirEdges2regionsOnOff = zeros(numEdgesDirectional,2);

edgeLIDs2regions = getRegionsForEdgeLIDs(setOfRegions_edgeLIDs,numEdges);
% numDirEdges = numEdges * 2;
for i=1:numRegions
    
    edgeLIDs_dir_region = c_edgeLIDsForRegions_cw{i};
    if(sum(edgeLIDs_dir_region>0)==0)
        continue
    end
    
    % get logical indices for edgeLID_dir
%     numDirEdge_region = numel(edgeLIDs_dir_region);
    % edgeLID_dir_sequence = 1:numDirEdges;
%     [~,edgeLIDs_dir_region_logical] = ...
%                 intersect(edgeLID_dir_sequence,edgeLIDs_dir_region);
    
    dirEdges2regionsOnOff(edgeLIDs_dir_region,1) = i; % =rID_on
    
    % get complementary edgeLIDs_dir for the set of edges
    edgeLID_dir_complementary = getComplementaryEdgeLID_dir...
                    (numEdges,edgeLIDs_dir_region);
    % for the complementary directions of these edges, region i is off
    dirEdges2regionsOnOff(edgeLID_dir_complementary,2) = i; % =rID_off
    % get rID_off
    % rID_off = getOffRID(edgeLIDs2regions,i);
    % dirEdges2regionsOnOff(edgeLIDs_dir_region,2) = rID_off; % rID_off  

end

function edgeLID_dir_complementary = getComplementaryEdgeLID_dir...
                    (numEdges,edgeLIDs_dir_region)
numRegionEdges = numel(edgeLIDs_dir_region);

edgeLID_dir_complementary = zeros(numRegionEdges,1);

addOffSetInds_logical = (edgeLIDs_dir_region<=numEdges);

edgeLID_dir_complementary(addOffSetInds_logical) = ...
        edgeLIDs_dir_region(addOffSetInds_logical) + numEdges;

subtractOffSetInds_logical = (edgeLIDs_dir_region>numEdges);
edgeLID_dir_complementary(subtractOffSetInds_logical) = ...
        edgeLIDs_dir_region(subtractOffSetInds_logical) - numEdges;
    
    

function rID_off = getOffRID(twoRegionEdgeLIDs,rID_on)
% col1
col1RID_on_logical = (twoRegionEdgeLIDs(:,1)==rID_on);
rID_off1 = twoRegionEdgeLIDs(col1RID_on_logical,2);

col2RID_on_logical = (twoRegionEdgeLIDs(:,2)==rID_on);
rID_off2 = twoRegionEdgeLIDs(col2RID_on_logical,1);

rID_off = [rID_off1; rID_off2];
