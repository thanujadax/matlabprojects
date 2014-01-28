function edges2nodes = getEdges2nodes_grid(adjacencyMat)

% get upper triangular matrix without diagonal (to avoid duplicates)
adjacencyMat = triu(adjacencyMat,1);  

% for each adjacency, define an edge
numEdges = sum(adjacencyMat>0);
edges2nodes = zeros(numEdges,2);

[edges2nodes(:,1),edges2nodes(:,2)] = find(adjacencyMat>0);