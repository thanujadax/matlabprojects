function [c_cells2regions,c_cellInternalEdgeIDs,c_cellIntNodeListInds] = getRegionsForOnCells(...
    faceAdj,onRegionIndList,offEdgeIDList,setOfRegions,wsIDsForRegions,ws,...
    offNodeIndList,edges2nodes,edgeListInds)

% Input:
%   faceAdj: face-adjacency graph for all the regions. The values are
%   edgeIDs: 
%   offEdgeIDList: 

% Outputs:

usedRegionsList = [];
numActiveRegions = numel(onRegionIndList);
k = 0; 

for i=1:numActiveRegions
    % check if region is already used
    if(sum(ismember(usedRegionsList,onRegionIndList(i)))==0)
        % if not used
        % regionList_i contains the regionIDs that are connected to each other
        regionList_i = [];
        internalEdgeList_i = [];
        [regionList_i,internalEdgeList_i] = getRegionList(...
            onRegionIndList(i),faceAdj,offEdgeIDList,regionList_i,internalEdgeList_i,...
            onRegionIndList,wsIDsForRegions,ws);
        usedRegionsList = [usedRegionsList regionList_i];
        k = k + 1;
        
        % regions
        c_cells2regions{k} = regionList_i;      
        % internal edges
        
        allEdgeIDsForRegions = setOfRegions(regionList_i,:);
        allEdgeIDsForRegions = allEdgeIDsForRegions(allEdgeIDsForRegions>0);
        internalEdgeList_i = intersect(offEdgeIDList,allEdgeIDsForRegions);
        
        c_cellInternalEdgeIDs{k} = internalEdgeList_i;
        
        % internal nodes
        edgeListInds_i = find(ismember(edgeListInds,internalEdgeList_i));
        nodeListIndsOfIntEdges_all = edges2nodes(edgeListInds_i,:);
        nodeListIndsOfIntEdges_all = nodeListIndsOfIntEdges_all(nodeListIndsOfIntEdges_all>0);
        nodeListInds_internal_i = intersect(offNodeIndList,nodeListIndsOfIntEdges_all);
        c_cellIntNodeListInds{k} = nodeListInds_internal_i; 
        
    end
end



function [thisRegionNeighborList,connectedEdgeIDs] = getRegionList(...
        thisRegionID,faceAdj,offEdgeIDList,thisRegionNeighborList,...
        connectedEdgeIDs,onRegionIndList,wsIDsForRegions,ws)

if(~(sum(ismember(thisRegionNeighborList,thisRegionID))))
thisRegionNeighborList = [thisRegionNeighborList thisRegionID];
end
% get all edges connected to this region
edgeIDsForThisRegion = faceAdj(thisRegionID,:); 
edgeIDsForThisRegion_nz = edgeIDsForThisRegion(edgeIDsForThisRegion>0);

% get active regions connected via off edges
offEdgeIDsForThisRegion = intersect(edgeIDsForThisRegion_nz,offEdgeIDList);
immediateNeighborRList = find(ismember(edgeIDsForThisRegion,offEdgeIDsForThisRegion));
activeNeighborRList = intersect(onRegionIndList,immediateNeighborRList);

if(~isempty(activeNeighborRList))
    newNeighbors = setdiff(activeNeighborRList,thisRegionNeighborList);
    if(~isempty(newNeighbors))
        if(size(newNeighbors,1)>1)
            newNeighbors = newNeighbors';
        end
        thisRegionNeighborList = [thisRegionNeighborList newNeighbors];
        newConnectedEdges = setdiff(offEdgeIDsForThisRegion,connectedEdgeIDs);
        if(~isempty(newConnectedEdges))
            if(size(newConnectedEdges,2)>1)
                newConnectedEdges = newConnectedEdges';
            end
            connectedEdgeIDs = [connectedEdgeIDs; newConnectedEdges];
        end
        % get the connected regions for the regions in the list as well
        % iterate until no new regions are added to the list of regions
        numNewNeighbors = numel(newNeighbors);
        for j=1:numNewNeighbors
            [thisRegionNeighborList,connectedEdgeIDs] = getRegionList(...
            newNeighbors(j),faceAdj,offEdgeIDList,thisRegionNeighborList,...
            connectedEdgeIDs,onRegionIndList,wsIDsForRegions,ws);
        end
    end
end
