function edgeUnary = getEdgeProbabilityFromMap(...
        membraneProbabilityImage,edgePixels,addMarginThickness,addMarginValue)
    
numEdges = size(edgePixels,1);
edgeUnary = zeros(numEdges,1);

% calculates edge unary probabilities from input membrane probability map

membraneProbabilityMap = double(imread(membraneProbabilityImage));
membraneProbabilityMap = membraneProbabilityMap./(max(max(membraneProbabilityMap)));
% add border to be compatible with filtered image
if(addMarginThickness>0)
    probabilityMapWithMargin = addThickBorder...
                (membraneProbabilityMap,addMarginThickness,addMarginValue);
else
    probabilityMapWithMargin = membraneProbabilityMap;
end

% get the edgeUnary to be the average edgeProbability (membrane
% probability)

for i=1:numEdges
    edgePix_i = edgePixels(i,:);
    edgePix_i = edgePix_i(edgePix_i>0);
    edgeUnary(i) = mean(probabilityMapWithMargin(edgePix_i));
end