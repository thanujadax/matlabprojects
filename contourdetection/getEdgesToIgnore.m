function edges2ignore = getEdgesToIgnore(edges2pixels,connectedJunctionIDs,...
            sizeR,sizeC)
% returns single pixel edges that only join already directly connected
% nodes

numEdges = size(edges2pixels,1);
numPixelsPerEdge = zeros(numEdges,1);
for i=1:numEdges
    numPixelsPerEdge(i) = sum(edges2pixels(i,:)>0);
end

singlePixelEdges = find(numPixelsPerEdge==1);
% check if all neighbors of any such single pixel edge are a clustered
% group of nodes
% get neighbors of the single pixel edges
numSinglePixEdges = numel(singlePixelEdges);
for i=1:numSinglePixEdges
    pixID = edges2pixels(singlePixelEdges(i));
    neighbors = getNeighbors(pixID,sizeR,sizeC);
    % see if the neighbors are in connectedJunctionIDs
    
end