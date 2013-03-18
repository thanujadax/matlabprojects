function nodeAngleCost = getNodeAngleCost(dTheta,sigma,C)
% calculate the cost at each node for each combination of edge pair. the
% array dTheta already gives the angle difference for each edge pair at
% each node in order

% C - parameter to scale the cost

% [numNodes,numAngles] = size(dTheta);
nodeAngleCost = gauss1d(dTheta,mean,sigma) .*C;
