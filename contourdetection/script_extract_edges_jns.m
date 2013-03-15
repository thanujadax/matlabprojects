% extracting edges and junctions from WS
% imIn = imread('stem_256x_t02_V.png');
imIn = imread('testMem4_V.png');
[sizeR,sizeC] = size(imIn);
ws = watershed(imIn);
%figure;imagesc(ws)

% % test input ws
% ws = imread('toyWS.png');
% [sizeR,sizeC] = size(ws);

% the edges (watershed boundaries) are labeled 0
% extract those
ind0 = find(ws==0);
% [r0,c0] = ind2sub([sizeR sizeC],ind0);
wsBoundaries = zeros(sizeR,sizeC);
wsBoundaries(ind0) = 1;
figure;imagesc(wsBoundaries);
title('watershed boundaries')

%% extracting junctions
% look at the 4 neighborhood of each pixel
fourNH = zeros(size(ws));
numBoundaryPixels = numel(ind0);  % watershed edge pixels
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
ind4J = find(fourNH>2);         % indices of junctions wrt to ws segmentation
% visualize junctions
wsJ = zeros(sizeR,sizeC);
wsJ(ind4J) = 1;
wsVis = zeros(sizeR,sizeC,3);
wsVis(:,:,3) = wsBoundaries;
wsVis(:,:,1) = wsJ;
figure;imagesc(wsVis);
title('Junctions from WS')

%% extracting edges connecting junctions
% assign unique labels for edges
% all the pixels in each edge should have the same label
wsEdges = wsBoundaries;
wsEdges(ind4J) = 0;         % setting junctions to zero. only edges are 1
pixList = find(wsEdges);    % edge pixels without junctions.
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
% TODO: assign colors from the OFR to each edge
edgePixColors = edgePixLabels;
edgePixColors(:,2) = mod(edgePixColors(:,2),100);
wsEdges2 = wsBoundaries;
wsEdges2(edgePixColors(:,1)) = edgePixColors(:,2);
figure;imagesc(wsEdges2);title('edges between junctions labeled separately')
%% extract edges with zero pixel length
connectedJunctionIDs = getClusteredJunctions(wsJ);
% connectedJunctionIDs contain the same ID for each node that is connected
% together with zero length edges
% nodeZeroEdges - store node - edge1,edge2 etc for these zero length edge
%% Build the adjacency matrix of the junction nodes
% for each node, get a list of edge IDs connected to it
[nodeEdges,nodeInds] = getNodeEdges(ind4J,edgePixLabels,connectedJunctionIDs,sizeR,sizeC);
[adjacencyMat,edges2nodes] = getAdjacencyMat(nodeEdges);
% visualize graph
% binary adjacency matrix
% [r,c] = ind2sub([sizeR sizeC],ind4J);
% k = 1:30;
% [B,XY] = bucky;
% gplot(B(k,k),XY(k,:),'-*')
% axis square

% edge to pixel correspondence
edges2pixels = getEdges2Pixels(edgePixLabels);

[r,c] = ind2sub([sizeR sizeC],nodeInds);
xy = [c r];
figure;gplot(adjacencyMat,xy,'-*');
set(gca,'YDir','reverse');
axis square

% n = size(adjacencyMat,1);
% k = 1:n;
% figure;gplot(adjacencyMat(k,k),xy(k,:),'-*')
% set(gca,'YDir','reverse');
% axis square

