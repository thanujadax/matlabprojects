function nodePixels = getNodeEdgePixel(nodeInd,edgePixelInds,sizeR,sizeC,maxNumPix)
% nodeInd contains the pixel index of a junction node
% edgePixelInds contains a set of pixel indices of an edge attached to the
% given node
% nodePixels (output) are the 3 edge pixel indices which are closest to the node
% pixel. If there are less than 3 such pixels, only those are returned.

if nargin < 5
    maxNumPix = 5;
end
[nodeR,nodeC] = ind2sub([sizeR sizeC],nodeInd);
[pixR,pixC] = ind2sub([sizeR sizeC],edgePixelInds);
numEdgePixels = numel(pixR);

if(numEdgePixels<maxNumPix)
    maxNumPix = numEdgePixels;
end
% select the 3 pixels which are closest to the junction node
distance = (pixR-nodeR).^2 + (pixC-nodeC).^2;
[sortedDists, sortedInd] = sort(distance);

sortedPixInds = edgePixelInds(sortedInd);

nodePixels = sortedPixInds(1:maxNumPix); 

