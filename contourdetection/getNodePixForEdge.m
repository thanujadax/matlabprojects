function nodePixForEdge = getNodePixForEdge(nodePixels,edgePixels,sizeR,sizeC)

% given a set of nodePixels (probably a cluster node) and a set of
% edgepixels, select to which nodePixel the edge is attached to. i.e. the
% closest

% for each edge pixel get the closest node pixel
% which distance is the least

nodePixForEdge = 0;

numEdgePixels = numel(edgePixels);
numNodePixels = numel(nodePixels);

distMatrix = zeros(numNodePixels,numEdgePixels);

for i=1:numNodePixels
    distances = getDistance(nodePixels(i),edgePixels,sizeR,sizeC);
    distMatrix(i,:) = distances;
end

[nodePixListInd,edgePixListInd] = find(distMatrix==min(distMatrix(:)));

if(numel(nodePixListInd)==1)
    nodePixForEdge = nodePixels(nodePixListInd);
else
    error('more than one node pixel found for given edge')
end
   


function distance = getDistance(nodePixInd,edgePixInds,sizeR,sizeC)
[refR,refC] = ind2sub([sizeR sizeC],nodePixInd);

[pixListR, pixListC] = ind2sub([sizeR sizeC],edgePixInds);

distance = (pixListR - refR).^2 + (pixListC - refC).^2;