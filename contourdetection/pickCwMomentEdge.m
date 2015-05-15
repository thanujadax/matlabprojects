function nextCwEdgeLId_inRegion = pickCwMomentEdge...
    (nextEdgeLIDs,prevEdgeLID,edges2pixels,regionID,ws)


% get the pixels around previous edge
% which of them belongs to this region?
% take the avg pixel location
% for each nextEdge, get the moment (angle) from the avg pix location
% retain the edge which gives a clockwise moment

% get pixels in region next to prevEdge
edgepixels = edges2pixels;
edgepixels(:,1) = [];

prevEdgePixels = edgepixels(preEdgeLID,:);
prevEdgePixels = prevEdgePixels(prevEdgePixels>0);

[sizeR,sizeC] = size(ws);
pixelsAroundPrevEdge = get8Neighbors(prevEdgePixels,sizeR,sizeC);
regionPixels = find(ws==regionID);
regionPixelsAroundEdge = intersect(regionPixels,pixelsAroundPrevEdge);

if(isempty(regionPixelsAroundEdge))
    error('no region pixels found!!')
end

medianRegionPixel = median(regionPixelsAroundEdge);