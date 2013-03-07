function adjacencyMat = getAdjacencyMat(nodeEdges)
% Input:
%   nodeEdges: gives the list of edgeIDs connected to each junction node
[numNodes, numEdgesPerNode] = size(nodeEdges);
adjacencyMat = zeros(numNodes);

numEdges = max(nodeEdges(:,2:numEdgesPerNode));

for i=1:numEdges
    % for each edge, find the two corresponding nodes at its ends
    [R,C] = find(nodeEdges(:,2:numEdgesPerNode));
    % R has the list indices of the junctions corresponding to edge i
    if(~isempty(R))
        % R not empty. assign to adjacencyMat
        nodeInd = nodeEdges(R,1); 
        adjacencyMat(nodeInd(1),nodeInd(2)) = i; % assign edgeId to the adjMat
    end    
end