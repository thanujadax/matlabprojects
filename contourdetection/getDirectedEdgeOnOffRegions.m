function dirEdges2regionsOnOff ...
        = getDirectedEdgeOnOffRegions...
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

% dirEdges2regionsOnOff - edgeListInd_dir (=rowID) | onRegion | offRegion  : dir N1->N2
%   regionID = 0 is for the image border.


% intermediate variables: (set as output if required. not need now) 
% c_edgeLIDsForRegions_dir_cw - set of cells each containing directional
%   edgeLIDs for the corresponding region
% setOfRegions_edgeLIDs - edgeLIDs for each region (undirected)
% edgeLIDs2nodes_directional - 

% for each edge, get the regions on either side.
% assing edgeLID for one region (on region)
% assign complementary edgeID for the other one (off region)
% for the two regions, get the directedEdgeLIDs on order.
% continue until there's no edge untraversed

edges2nodes_complements = edges2nodes;
edges2nodes_complements(:,1) = edges2nodes(:,2);
edges2nodes_complements(:,2) = edges2nodes(:,1);
edgeLIDs2nodes_directional = [edges2nodes; edges2nodes_complements];

dirEdges2regionsOnOff = zeros(size(edgeLIDs2nodes_directional,1),2);

untraversedEdgeLIDs = edgeListIndsAll;

while(numel(untraversedEdgeLIDs)>0)
    edgeID_i = untraversedEdgeLIDs(end);
    untraversedEdgeLIDs(end) = []; % removed the current edge from the untr
    edgeLID_i = find(edgeListIndsAll==edgeID_i);
    complementaryEdgeLID_i = edgeLID_i + numel(edgeListIndsAll);
    % get the regions
    [regionIDs,~] = find(setOfRegions==edgeID_i);
    % if there's only one region, the other region is the border rID = 0;
    if(isempty(regionIDs))
        error('no regions found for edge!')
    elseif(numel(regionIDs)==1)
        regionIDs(end+1) = 0;
    elseif(numel(regionIDs)>2)
        error('more than 2 regions found for edge!!')
    end
    
    % on side for edge
    dirEdges2regionsOnOff(edgeLID_i,1) = regionIDs(1);
    % get all the other directedEdgeLIDs for the same region in the same direction
    % set current region as on region for these edges
    dirEdgeLIDsForRegion = getDirEdgeLIDsForRegion...
            (regionIDs(1),edgeLID_i,setOfRegions,edgeLIDs2nodes_directional,...
            edgeListIndsAll);
    % get the other region for each complementary edge and set them as off
    % regions
    % remove these regions from the untraversed list of edges
    
    % off side for edge
    dirEdges2regionsOnOff(complementaryEdgeLID_i,2) = regionIDs(2);
    
    
    
end

