function [c_cells2regions,c_cellInternalEdgeIDs] = getRegionsForOnCells(faceAdj,...
            onRegionIndList,offEdgeIDList,setOfRegions,wsIDsForRegions,ws)

% Input:
%   faceAdj: face-adjacency graph for all the regions. The values are
%   edgeIDs: 
%   offEdgeIDList: 

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
        c_cells2regions{k} = regionList_i;
        c_cellInternalEdgeIDs{k} = internalEdgeList_i;
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
        thisRegionNeighborList = [thisRegionNeighborList newNeighbors];
        newConnectedEdges = setdiff(offEdgeIDsForThisRegion,connectedEdgeIDs);
        if(~isempty(newConnectedEdges))
            connectedEdgeIDs = [connectedEdgeIDs; newConnectedEdges'];
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
