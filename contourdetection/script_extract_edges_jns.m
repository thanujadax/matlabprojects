% extracting edges and junctions from WS
imIn = imread('stem_256x_t02_V.png');
[sizeR,sizeC] = size(imIn);
ws = watershed(imIn);
%figure;imagesc(ws)

% the edges (watershed boundaries) are labeled 0
% extract those
ind0 = find(ws==0);
% [r0,c0] = ind2sub([sizeR sizeC],ind0);
wsBoundaries = zeros(sizeR,sizeC);
wsBoundaries(ind0) = 1;
figure;imshow(wsBoundaries);
title('watershed boundaries')

%% extracting junctions
% look at the 4 neighborhood of each pixel
fourNH = zeros(size(imIn));
numBoundaryPixels = numel(ind0);
for i=1:numBoundaryPixels
    % calculate n.o. 4 neighbors
    ind = ind0(i);
    [r c] = ind2sub([sizeR sizeC],ind);
    nh = zeros(1,4);
    
    if((r-1)>0)
        nh(1) = wsBoundaries(r-1,c);
    end
    if((r+1)<=sizeR)
        nh(2) = wsBoundaries(r+1,c);
    end
    if((c-1)>0)
        nh(3) = wsBoundaries(r,c-1);
    end
    if((c+1)<=sizeC)
        nh(4) = wsBoundaries(r,c+1);
    end
    
    fourNH(ind) = sum(nh);
end
% get the pixels which are having a 4NH > 2
ind4J = find(fourNH>2);         % indices of junctions
% visualize junctions
wsJ = zeros(sizeR,sizeC);
wsJ(ind4J) = 1;
wsVis = zeros(sizeR,sizeC,3);
wsVis(:,:,3) = wsBoundaries;
wsVis(:,:,1) = wsJ;
figure;imshow(wsVis);
title('Junctions from WS')

%% extracting edges connecting junctions
% assign unique labels for edges
% all the pixels in each edge should have the same label
wsEdges = wsBoundaries;
wsEdges(ind4J) = 0;         % setting junctions to zero. only edges are 1
pixList = find(wsEdges);
edgePixLabels = zeros(numel(pixList),2); % 2nd column stores the edge label
edgePixLabels(:,1) = pixList;
currentLabel = 0;
for i=1:numel(pixList)
    if(edgePixLabels(i,2)==0)
        % assign label
        currentLabel = currentLabel + 1;
        edgePixLabels = labelPixelNeighbors(edgePixLabels,pixList(i),currentLabel,...
                sizeR,sizeC);
    else
        continue
    end
end
% assign random colors to edges
edgePixColors = edgePixLabels;
edgePixColors(:,2) = mod(edgePixColors(:,2),100);
wsEdges2 = wsBoundaries;
wsEdges2(edgePixColors(:,1)) = edgePixColors(:,2);
figure;imagesc(wsEdges2);title('edges between junctions labeled separately')

%% Build the adjacency matrix of the junction nodes
% for each node, get a list of edge IDs connected to it
nodeEdges = getNodeEdges(ind4J,edgePixLabels,sizeR,sizeC);
adjacencyMat = getAdjacencyMat(nodeEdges);




