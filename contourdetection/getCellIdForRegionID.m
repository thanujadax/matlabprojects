function cellID = getCellIdForRegionID(regionID, c_cells2regions)
% get cell given region
cellID = 0;
numCells = numel(c_cells2regions);
for i=1:numCells
    setOfRegions_i = c_cells2regions{i};
    if(sum(ismember(setOfRegions_i,regionID)))
        % region belongs to cell i
        cellID = i;
        % break
    end
end
