function c_cells2regions = getRegionsForCells(faceAdj,offEdgeIDList)

usedRegionList = [];
numRegions = size(faceAdj,1);

for i=1:numRegions
    % regionList_i contains the regionIDs that are connected to each other
    [regionList_i,usedRegionList] = getRegionList(faceAdj,offEdgeIDList); % recursive function
    
end



function [regionList,usedRegionList] = getRegionList(faceAdj,offEdgeIDList)
numRegions = size(faceAdj,1);
% for each region

% get all edges connected to it
% get edges turned off
% get regions connected via off edges

