function c_cells2regions = getRegionsForCells(faceAdj,offEdgeIDList)

usedRegionsList = [];
numRegions = size(faceAdj,1);
k = 0;
for i=1:numRegions
    % check if region is already used
    if(sum(ismember(usedRegionsList,i))==0)
        % if not used
        % regionList_i contains the regionIDs that are connected to each other
        clear regionList_i
        [regionList_i,usedRegionsList] = getRegionList(i,faceAdj,offEdgeIDList,usedRegionsList);
        k = k + 1;
        c_cells2regions{k} = regionList_i;
    end
end



function [regionList,usedRegionsList] = getRegionList(thisRegionID,faceAdj,offEdgeIDList,...
                usedRegionsList)
numRegions = size(faceAdj,1);
allRegionsList = 1:numRegions;

% get all edges connected to this region
edgeIDsForThisRegion = faceAdj(thisRegionID,:); 
edgeIDsForThisRegion_nz = edgeIDsForThisRegion(edgeIDsForThisRegion>0);

% get regions connected via off edges
offEdgeIDsForThisRegion = intersect(edgeIDsForThisRegion_nz,offEdgeIDList);
regionList=find(ismember(edgeIDsForThisRegion,offEdgeIDsForThisRegion));

usedRegionsList = [usedRegionsList thisRegionID regionList];
