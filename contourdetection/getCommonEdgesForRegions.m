function edgeLIDsCommonToMergedRegions = getCommonEdgesForRegions(...
    listOfMergedRegions,ws,edges2pixels)

% For a given set of regions, consider each pair and figure out which edges
% are interfacing both regions

% Inputs:
%   listOfMergedRegions - contains old regionIDs
%   ws - original ws
%   edges2pixels

% Outputs:
%   edgeListInds of the common edges for the region list given

edgeLIDsCommonToMergedRegions = [];

cell_regionEdgeLIDs = {};
numMergedRegions = numel(listOfMergedRegions);

for i=1:numMergedRegions
    % get edges for each region
    cell_regionEdgeLIDs{i} = getRegionEdges(listOfMergedRegions(i),ws,edges2pixels);
end

% see if there are common edges for each pair of merged regions
% get all pair combinations of regions
V = 1:numMergedRegions;
if(numMergedRegions>1)
    % if there's only one region, we don't have to do anything since there
    % won't be any additional edges which are potentially dividing multiple
    % regions that were merged afterwards
    
    regionCombinationInds = nchoosek(V,2);
    numCombinations = size(regionCombinationInds,1);
    % for each combination (pair of regions) check if there is a common edge

    for i=1:numCombinations
        regionLID_a = regionCombinationInds(i,1);
        edgeLIDs_a = cell_regionEdgeLIDs{regionLID_a};

        regionLID_b = regionCombinationInds(i,2);
        edgeLIDs_b = cell_regionEdgeLIDs{regionLID_b};    

        commonEdgeLIDs = intersect(edgeLIDs_a,edgeLIDs_b);

        edgeLIDsCommonToMergedRegions = [edgeLIDsCommonToMergedRegions; commonEdgeLIDs];
    end
end


function regionEdgeLIDs = getRegionEdges(regionID,ws,edges2pixels)

intPix = (ws==regionID);
regionEdgePix = getEdgePixForWsFace(intPix,ws);
% get the edgeLIDs for the pixels
[sizeR,sizeC] = size(ws);
regionEdgeLIDs = getEdgeIDsFromPixels(regionEdgePix,edges2pixels);