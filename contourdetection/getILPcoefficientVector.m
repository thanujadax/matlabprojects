function f = getILPcoefficientVector(edgePriors,j3NodeAngleCost,j4NodeAngleCost)
numEdges = size(edgePriors,1);
numJ3 = size(j3NodeAngleCost,1);
numJ4 = size(j4NodeAngleCost,1);

numElements = numEdges*2 + numJ3*4 + numJ4*7;
f = zeros(numElements,1);
% order of elements
%[{edgeInactive},{edgeActive},{J3inactive},{J3Active_3},{J4inactive},{J4Active_6}]

% edge variables
j=1;
for i=1:2:2*numEdges
    f(i) = edgePriors(j);       % inactivation cost
    f(i+1) = -edgePriors(j);  % activation cost for the same edge
    j = j+1;
end

% junction variables - J3
maxJ3cost = max(j3NodeAngleCost,[],2);
% j3NodeAngleCost = -j3NodeAngleCost;
j3NodeAngleCost = [maxJ3cost j3NodeAngleCost]; % first column: inactivation cost

numJ3Coeff = numel(j3NodeAngleCost);
f((2*numEdges+1):(2*numEdges+numJ3Coeff)) = j3NodeAngleCost(1:numJ3Coeff);

% junction variables - J4
maxJ4cost = max(j4NodeAngleCost,[],2);
j4NodeAngleCost = -j4NodeAngleCost;
j4NodeAngleCost = [maxJ4cost j4NodeAngleCost]; % first column: inactivation cost

numJ4Coeff = numel(j4NodeAngleCost);
f((2*numEdges+numJ3Coeff+1):(2*numEdges+numJ3Coeff+numJ4Coeff)) = ...
                j4NodeAngleCost(1:numJ4Coeff);
