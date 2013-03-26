function nodeAngleCost = getNodeAngleCostExp(dTheta,midpoint,rateConst,C)
% calculate the cost at each node for each combination of edge pair. the
% array dTheta already gives the angle difference for each edge pair at
% each node in order

% C - parameter to scale the cost

% [numNodes,numAngles] = size(dTheta);
numJ = numel(dTheta);
[rows,cols] = size(dTheta);
nodeAngleCost = zeros(rows,cols);
for i=1:rows
    for j=1:cols
        t = dTheta(i,j);
        if(t<0)
            % invalid angle difference
            nodeAngleCost(i,j) = 1;
        else
            % nodeAngleCost = gauss1d(dTheta,midpoint,rateConst) .*C;
            if(t<180)
                t = midpoint*2 - t;
            end        
            x = t - midpoint;
            exponent = -1*rateConst*x;
            nodeAngleCost(i,j) = exp(exponent)*C;    
        end
    end
end