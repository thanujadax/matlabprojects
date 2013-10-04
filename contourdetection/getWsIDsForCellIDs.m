function wsIDs = getWsIDsForCellIDs(ws,setOfCells,edges2pixels,nodeInds,...
            edges2nodes,edgeIDsAll)
% Output:
%   wsIDs - each element gives the ID in WS, in the order of the list of
%   cells given in setOfCells

% PROBLEM: haven't included clustered junctions

numCells = size(setOfCells,1);
[sizeR,sizeC] = size(ws);
wsIDs = zeros(numCells,1);

for i=1:numCells
    % get internal pixels
    edgeSet_cell = setOfCells(i,:);
    edgeSet_cell = edgeSet_cell(edgeSet_cell>0);
    boundaryPixels = getBoundaryPixelsForCell(edgeSet_cell,edges2pixels,...
    nodeInds,edges2nodes,edgeIDsAll);
    [internalx,internaly] = getInternelPixelsFromBoundary(boundaryPixels,sizeR,sizeC);
    % get the most popular label for the internal pixels
    internalPixels = sub2ind([sizeR sizeC],internaly,internalx);
    internalPixWsLabels = ws(internalPixels);
    internalPixWsLabels = internalPixWsLabels(internalPixWsLabels>0);
    wsIDs(i) = mode(double(internalPixWsLabels));   
end