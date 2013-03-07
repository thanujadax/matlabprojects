function nodeEdges = getNodeEdges(nodeInd,edgePixLabels,sizeR,sizeC)
% Inputs:
%   nodeInd - array of junction indices
%   edgePixLabels - N-by-2 array of edge labels for each pixel, given by
%   the index wrt the original image (watershed)

% Output:
%   nodeEdges - array with edge labels corresponding to each junction. each
%   row -> jn

% for each node, get the neighbors
% get the edgeID of the neighbors
nodeEdges = zeros(numel(nodeInd),4);

for i=1:numel(nodeInd)
    neighborInd = getNeighbors(nodeInd(i),sizeR,sizeC);
    numNeighbors = numel(neighborInd);
    k = 0;
    for j=1:numNeighbors
        neighborListInd = find(edgePixLabels(:,1)==neighborInd(j));
        if(~isempty(neighborListInd))
            % there's an edge for this neighbor
            k = k + 1;
            nodeEdges(i,k) = edgePixLabels(neighborListInd,2);
        end
    end
end