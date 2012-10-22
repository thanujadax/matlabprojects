function adj = getAdjMatrix(NN1,NN2,bb)
% returns upper triangular adjacency matrix

numNodesHorizontal = NN - bb + 1;
numNodesVertical = NN - bb + 1;
numNodes = numNodesHorizontal * numNodesVertical;

adj = zeros(numNodes);

% the adjacency matrix is build up from left to right and top to bottom
% as an upper triangular matrix
for i = 1 : numNodes
    j = i + 1;    
    if (mod(i,numNodesHorizontal) > 0)  % excludes the last column
        adj(i,j) = 1;           % neighbor to the right
    end
    j = i + numNodesHorizontal
    if i <= (numNodes - numNodesHorizontal) % excludes the last row
        adj(i,j) = 1;           % neighbor to the bottom
    end
end


    

