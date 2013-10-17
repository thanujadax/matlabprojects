function regionID = getRegionIdForEdgeID(edgeID, c_setOfRegions)
% get regionID given edgeID
regionID = 0;
numRegions = numel(c_setOfRegions);
for i=1:numRegions
    setOfRegions_i = c_setOfRegions{i};
    if(sum(ismember(setOfRegions_i,edgeID)))
        % region belongs to cell i
        regionID = i;
        % break
    end
end