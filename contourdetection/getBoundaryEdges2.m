function boundaryEdges = getBoundaryEdges2(wsgraph,marginSize,edges2pixels,...
    edges2nodes,nodeEdges,edgeIDs)

% from each boundary, along perpendicular lines get all the pixels nearest
% to the edge i.e. along each line given by a set of pixels.
% Remove the nodeInds from this list of pixels. Then we are left with the
% edge pixels. Get the boundaryEdgeIDs from these

[sizeR, sizeC] = size(wsgraph); % size of the image
% top edge pixels
topEdgePixels = zeros(sizeC,1);
for i=1:sizeC
    for j=1:marginSize
        if(wsgraph(j,i)>0)
            topEdgePixels(i) = sub2ind([sizeR sizeC],j,i);
            break
        end
    end
end

% bottom edge pixels
bottomEdgePixels = zeros(sizeC,1);
for i=1:sizeC
    for j=sizeR:-1:(sizeR-marginSize)
        if(wsgraph(j,i)>0)
            bottomEdgePixels(i) = sub2ind([sizeR sizeC],j,i);
            break
        end
    end
end

% left edge pixels
leftEdgePixels = zeros(sizeR,1);
for i=1:sizeR
    for j=1:marginSize
        if(wsgraph(i,j)>0)
            leftEdgePixels(i) = sub2ind([sizeR sizeC],i,j);
            break
        end
    end
end

% right edge pixels
rightEdgePixels = zeros(sizeR,1);
for i=1:sizeR
    for j=sizeC:-1:(sizeC-marginSize)
        if(wsgraph(i,j)>0)
            rightEdgePixels(i) = sub2ind([sizeR sizeC],i,j);
            break
        end
    end
end

allBoundaryEdgePixels = [topEdgePixels; bottomEdgePixels;...
                leftEdgePixels; rightEdgePixels];
allBoundaryEdgePixels = unique(allBoundaryEdgePixels);            

% remove the nodeInds
nodeInds = nodeEdges(:,1);
boundaryEdgePixels = setdiff(allBoundaryEdgePixels,nodeInds);

% for each boundaryEdgePixel, get the corresponding edgeID
numBoundaryEdgePixels = numel(boundaryEdgePixels);
boundaryEdges = zeros(numBoundaryEdgePixels,1);
for i=1:numBoundaryEdgePixels
    clear edgeListInd
    [edgeListInd,~] = find(edges2pixels==boundaryEdgePixels(i));
    boundaryEdges(i) = edges2pixels(edgeListInd,1);
end

boundaryEdges = unique(boundaryEdges);