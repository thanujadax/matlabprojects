function [c_cells2regions,c_cellInternalEdgeIDs] = getRegionsForCells(faceAdj,offEdgeIDList)

% Input:
%   faceAdj: face-adjacency graph for all the regions. The values are
%   edgeIDs: 
%   offEdgeIDList: 

usedRegionsList = [];
numRegions = size(faceAdj,1);
k = 0;
for i=1:numRegions
    % check if region is already used
    if(sum(ismember(usedRegionsList,i))==0)
        % if not used
        % regionList_i contains the regionIDs that are connected to each other
        regionList_i = [];
        internalEdgeList_i = [];
        [regionList_i,internalEdgeList_i] = getRegionList(...
            i,faceAdj,offEdgeIDList,regionList_i,internalEdgeList_i);
        usedRegionsList = [usedRegionsList regionList_i];
        k = k + 1;
        c_cells2regions{k} = regionList_i;
        c_cellInternalEdgeIDs{k} = internalEdgeList_i;
    end
end



function [thisRegionNeighborList,connectedEdgeIDs] = getRegionList(...
        thisRegionID,faceAdj,offEdgeIDList,thisRegionNeighborList,...
        connectedEdgeIDs)
numRegions = size(faceAdj,1);

% get all edges connected to this region
edgeIDsForThisRegion = faceAdj(thisRegionID,:); 
edgeIDsForThisRegion_nz = edgeIDsForThisRegion(edgeIDsForThisRegion>0);

% get regions connected via off edges
offEdgeIDsForThisRegion = intersect(edgeIDsForThisRegion_nz,offEdgeIDList);
immediateNeighborList = find(ismember(edgeIDsForThisRegion,offEdgeIDsForThisRegion));

if(~isempty(immediateNeighborList))
    newNeighbors = setdiff(immediateNeighborList,thisRegionNeighborList);
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
            connectedEdgeIDs);
        end
    end
end
