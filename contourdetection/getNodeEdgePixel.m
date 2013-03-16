function nodePixel = getNodeEdgePixel(nodeInd,edgePixelInds,sizeR,sizeC)
% nodeInd contains the pixel index of a junction node
% edgePixelInds contains a set of pixel indices of an edge attached to the
% given node
% nodePixel (output) is the edge pixel index which is closest to the node
% pixel

[nodeR,nodeC] = ind2sub([sizeR sizeC],nodeInd);
[pixR,pixC] = ind2sub([sizeR sizeC],edgePixelInds);

distance = (pixR-nodeR).^2 + (pixC-nodeC).^2;
minDistance = min(distance);
nodePixListInd = find(distance==minDistance);
nodePixel = edgePixelInds(nodePixListInd);

