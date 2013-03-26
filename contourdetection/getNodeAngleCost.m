function nodeAngleCost = getNodeAngleCost(dTheta,midpoint,sigma,C)
% calculate the cost at each node for each combination of edge pair. the
% array dTheta already gives the angle difference for each edge pair at
% each node in order

% C - parameter to scale the cost

% [numNodes,numAngles] = size(dTheta);
if(dTheta<0)
    % invalid angle difference
    nodeAngleCost = 1;
else
    nodeAngleCost = gauss1d(dTheta,midpoint,sigma) .*C;
end
