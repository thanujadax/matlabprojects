function nodeAngleCost = getNodeAngleCost_directional(theta,alpha,C)
% calculate the cost for each active configuration of a node given the edge
% angles according to the OFR (theta) and the edge position of the graph
% relative the current node (anlpha)

% the number of columns of each theta and alpha corresponds to the number
% of edges connected to a node. this gives rise to a combination of the
% pairs of edges we can pick - nchoosek(n,2). in the output: nodeAngleCost,
% we have a number of columns equal to the number of such combinations

[numNodes,numEdgesPerNode] = size(theta);
numCombinations = nchoosek(numEdgesPerNode,2);
edgeIDvect = 1:numNodes;
combinations = nchoosek(edgeIDvect,2);

% get the outwardness score for all the edge for all the nodes
outwardnessScores = getOutwardness(theta,alpha,C);

% calculate the cost for each pair of edge combinations - nodeAngleCost
% = multiplication of the outwardness scores of the two edges
nodeAngleCost = zeros(numNodes,numCombinations);
for i=1:numNodes
   for j=1:numCombinations
        edge1LInd = combinations(j,1);
        edge2LInd = combinations(j,2);
        nodeAngleCost(i,j) = outwardnessScores(i,edge1LInd) *...
                            outwardnessScores(i,edge2LInd);  
        
   end
end
