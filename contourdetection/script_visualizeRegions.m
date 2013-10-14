intPixels = [];
[sizeR,sizeC]=size(ws);
thisRegionNeighborList = onRegionIndList;
numCells = numel(thisRegionNeighborList);
for i=1:numCells
    wsID_i = wsIDsForRegions(thisRegionNeighborList(i));
    if(wsID_i==0)
        disp('problem with wsid check')
    else

        cellPixInds_i = getInternalPixForCell(ws,wsID_i);
        intPixels = [intPixels; cellPixInds_i];
    end
end

regionSet = zeros(sizeR,sizeC);
regionSet(intPixels) = 1;
figure;imshow(regionSet)