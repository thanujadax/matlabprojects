function f = getILPcoefficientVector(edgePriors,j3NodeAngleCost,j4NodeAngleCost)
numEdges = size(edgePriors,1);
numJ3 = size(j3NodeAngleCost,1);
numJ4 = size(j4NodeAngleCost,1);

numElements = numEdges*2 + numJ3*4 + numJ4*7;
f = zeros(numElements,1);
% order of elements
%[{edgeInactive},{edgeActive},{J3inactive},{J3Active_3},{J4inactive},{J4Active_6}]

% edge variables
f(1:numEdges) = edgePriors;
f((numEdges+1):(2*numEdges)) = -edgePriors;

% junction variables - J3
maxJ3cost = max(j3NodeAngleCost,[],2);
j3NodeAngleCost = -j3NodeAngleCost;
j3NodeAngleCost = [maxJ3cost j3NodeAngleCost];
for i=(2*numEdges+1):(2*numEdges+numJ3)
    
end