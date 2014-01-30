function [edges2nodes,nodeEdges] = getEdges2nodes_grid(adjacencyMat)

% works for any type of grid

% get upper triangular matrix without diagonal (to avoid duplicates)
adjacencyMat_u = triu(adjacencyMat,1);  

% for each adjacency, define an edge
numEdges = sum(sum(adjacencyMat_u>0));
edges2nodes = zeros(numEdges,2);

[edges2nodes(:,1),edges2nodes(:,2)] = find(adjacencyMat_u>0);

% nodeEdges
numNodes = size(adjacencyMat_u,1);

nonZeros = adjacencyMat_u>0;
numNeighborsPerNode = sum(nonZeros);
maxEdgesPerNode = max(numNeighborsPerNode);

nodeEdges = zeros(numNodes,maxEdgesPerNode);
for i=1:numNodes
    nodeEdges_i = adjacencyMat(i,:);
    nodeEdges_i = nodeEdges_i(nodeEdges_i>0);
    numNodeEdges_i = numel(nodeEdges_i);
    if(numNodeEdges_i>0)
        nodeEdges(i,1:numNodeEdges_i) = nodeEdges_i;
    end
end