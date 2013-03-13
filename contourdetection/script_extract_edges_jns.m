% extracting edges and junctions from WS
% imIn = imread('stem_256x_t02_V.png');
imIn = imread('testMem4_V.png');
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
figure;imshow(wsVis);
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
% junction nodes (pixels) which are next to each other
% wsJ contains 1 for junction nodes, others zero.

eightNH_J = zeros(sizeR,sizeC);
numJunctionPixels = numel(ind4J);  % watershed junction pixels (3 and 4 Js)
for i=1:numJunctionPixels
    % calculate n.o. 4 neighbors
    ind = ind4J(i);
    [r c] = ind2sub([sizeR sizeC],ind);
    nh = zeros(1,8);
    % N1
    if((r-1)>0)
        nh(1) = wsBoundaries(r-1,c);
    end
    % N2
    if((r+1)<=sizeR)
        nh(2) = wsBoundaries(r+1,c);
    end
    % N3
    if((c-1)>0)
        nh(3) = wsBoundaries(r,c-1);
    end
    % N4
    if((c+1)<=sizeC)
        nh(4) = wsBoundaries(r,c+1);
    end
    % N5
    if((r-1)>0 && (c-1)>0)
        nh(5) = wsBoundaries(r-1,c-1);
    end
    % N6
    if((r-1)>0 && (c+1)<=sizeC)
        nh(6) = wsBoundaries(r-1,c+1);
    end
    % N7
    if((r+1)<=sizeR && (c-1)>0)
        nh(7) = wsBoundaries(r+1,c-1);
    end
    % N8
    if((r+1)<=sizeR && (c+1)<=sizeC)
        nh(8) = wsBoundaries(r+1,c+1);
    end
    eightNH_J(ind) = sum(nh);
end

indJClusterPixels = find(eightNH_J>0);  % gets junction pixels which has neighbors
% make a N-by-2 array of such junction pixels having immediate neighboring
% junctions
numJunctionClustPix = numel(indJClusterPixels);
connectedJunctionIDs = zeros(numJunctionClustPix);
connectedJunctionIDs(:,1) = indJClusterPixels;
junctionID = 0;
for i=1:numJunctionClustPix
    % look for the neighbors and give them the same label
    jLabelCurrent = connectedJunctionIDs(i,2);
    if(jLabelCurrent==0)
        % assign new label to this junction and all its neighbors
        junctionID = junctionID + 1;
        connectedJunctionIDs(i,2) = junctionID;
        % now deal with its neighbors
        ind = connectedJunctionIDs(i,1);
        [r c] = ind2sub([sizeR sizeC],ind);
        % N1
        if((r-1)>0)
            isJunction = eightNH_J((r-1),c);
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],(r-1),c);
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
        % N2
        if((r+1)<=sizeR)
            isJunction = eightNH_J((r+1),c);
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],(r+1),c);
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
        % N3
        if((c-1)>0)
            isJunction = eightNH_J(r,(c-1));
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],r,(c-1));
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
        % N4
        if((c+1)<=sizeC)
            isJunction = eightNH_J(r,(c+1));
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],r,(c+1));
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
        % N5
        if((r-1)>0 && (c-1)>0)
            isJunction = eightNH_J((r-1),(c-1));
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],(r-1),(c-1));
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
        % N6
        if((r-1)>0 && (c+1)<=sizeC)
            isJunction = eightNH_J((r-1),(c+1));
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],(r-1),(c+1));
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
        % N7
        if((r+1)<=sizeR && (c-1)>0)
            isJunction = eightNH_J((r+1),(c-1));
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],(r+1),(c-1));
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
        % N8
        if((r+1)<=sizeR && (c+1)<=sizeC)
            isJunction = eightNH_J((r+1),(c+1));
            if(isJunction)
                indNeigh = sub2ind([sizeR sizeC],(r+1),(c+1));
                listIndNeigh = find(connectedJunctionIDs(:,1));
                connectedJunctionIDs(listIndNeigh,2) = junctionID;
            end
        end
    end
end
% connectedJunctionIDs contain the same ID for each node that is connected
% together with zero length edges
% nodeZeroEdges - store node - edge1,edge2 etc for these zero length edges

%% Build the adjacency matrix of the junction nodes
% for each node, get a list of edge IDs connected to it
nodeEdges = getNodeEdges(ind4J,edgePixLabels,connectedJunctionIDs,sizeR,sizeC);
adjacencyMat = getAdjacencyMat(nodeEdges);
% visualize graph
% binary adjacency matrix
[r,c] = ind2sub([sizeR sizeC],ind4J);
xy = [c r];
figure;gplot(adjacencyMat,xy);
set(gca,'YDir','reverse');


