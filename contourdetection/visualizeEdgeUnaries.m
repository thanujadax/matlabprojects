function edgeUnaryMat = visualizeEdgeUnaries(edgepixels,edgePriors,sizeR,sizeC)
% Takes in the calculated edgeUnary values and assigns them the relevent
% edge pixels to create a monochromatic visualization of the edge unary
% intensities based on the abs OFR (or any other criteria used to calculate
% edgePriors

% normalize edgePriors if not normalized
maxUnary = max(edgePriors);
if(maxUnary>1)
    % normalize
    edgePriors = edgePriors./maxUnary;
end
numEdges = size(edgepixels,1);
edgeUnaryMat = zeros(sizeR,sizeC);
for i=1:numEdges
    edgepix_i = edgepixels(i,:);
    edgepix_i = edgepix_i(edgepix_i>0);
    edgeUnaryMat(edgepix_i) = edgePriors(i);
end