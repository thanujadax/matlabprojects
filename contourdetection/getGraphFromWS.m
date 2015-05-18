function [adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs,...
    selfEdgePixelSet,ws,ws_original,removedWsIDs,newRemovedEdgeLIDs]...
    = getGraphFromWS(ws,hsvOutput,displayImg)

% Outputs:
% nodeEdges: contains the set of edgeIDs for each nodePixInd

saveMatrices = 0;  % to save some of the generated matrices

% % test input ws
% ws = imread('toyWS.png');
[sizeR,sizeC] = size(ws);

% the edges (watershed boundaries) are labeled 0
% extract those
ind0 = find(ws==0);
% [r0,c0] = ind2sub([sizeR sizeC],ind0);
wsBoundaries = zeros(sizeR,sizeC);
wsBoundaries(ind0) = 1;
% figure;imshow(wsBoundaries);
% title('watershed boundaries')

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
% figure;imshow(wsVis);
% title('Junctions from WS')
% ws edges with OFR color code
hsvOutput_V = hsvOutput(:,:,3);
edgepix = zeros(sizeR,sizeC);
edgepix(wsBoundaries>0) = hsvOutput_V(wsBoundaries>0);
edgepix(wsJ>0) = 1;
hsvOutput(:,:,3) = edgepix;
if(saveMatrices)
    save('edgepix.mat','edgepix');
end
hsvImg = cat(3,hsvOutput(:,:,1),hsvOutput(:,:,2),hsvOutput(:,:,3));
RGBimg = hsv2rgb(hsvImg);
if(displayImg)
    figure;imshow(RGBimg);
end
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
%% visualization
% assign random colors to edges
% TODO: assign colors from the OFR to each edge
edgePixColors = edgePixLabels;
edgePixColors(:,2) = mod(edgePixColors(:,2),100);
wsEdges2 = wsBoundaries;
wsEdges2(edgePixColors(:,1)) = edgePixColors(:,2);
% figure;imshow(wsEdges2);title('edges between junctions labeled separately')
%% extract edges with zero pixel length
% connectedJunctionIDs = getClusteredJunctions(wsJ);
[connectedJunctionIDs,psuedoEdges2nodes] = getClusterNodesAndPsEdges(wsJ);
% connectedJunctionIDs contain the same ID for each node that is connected
% together with zero length edges, in the 4 neighborhood.
% pEdges2nodes - each row will give 2 nodes connected by a zero length edge

%% Build the adjacency matrix of the junction nodes
% edge to pixel correspondence
edges2pixels = getEdges2Pixels(edgePixLabels);
% edges2ignore = getEdgesToIgnore(edges2pixels,connectedJunctionIDs,sizeR,sizeC);
% for each node, get a list of edge IDs connected to it

numPsuedoEdges = size(pEdges2Nodes,1);
maxEdgeID = size(edges2pixels,1);
psuedoEdgeIDs = (maxEdgeID+1) : (maxEdgeID+numPsuedoEdges);

[nodeEdges,nodeInds] = getNodeEdges(ind4J,edgePixLabels,connectedJunctionIDs,sizeR,sizeC,...
            psuedoEdgeIDs,psuedoEdges2nodes);

[adjacencyMat,edges2nodes,selfEdgeIDs,~] = getAdjacencyMat(nodeEdges);

% calculate new ws by merging those ws regions that were initially separated
ws_original = ws;
[ws,removedWsIDs, newRemovedEdgeLIDs] = getCorrectedWSregions(ws,selfEdgeIDs,edges2pixels,displayImg);

% initialize output
selfEdgePixelSet = zeros(numel(selfEdgeIDs),1);

if(selfEdgeIDs(1)~=0)
    % remove selfEdges from nodeEdges, edges2nodes and edges2pixels
    % edges2nodes
    edges2nodes = edges2nodes((edges2nodes(:,1)~=0),:);

    [nodeEdgeRows,nodeEdgeCols] = size(nodeEdges);
    numSelfEdges = numel(selfEdgeIDs);
    for i=1:numSelfEdges
        % nodeEdges
        [rx,cx] = find(nodeEdges(:,2:nodeEdgeCols)==selfEdgeIDs(i));
        cx = cx + 1;
        nodeEdges(rx,cx)=0;    
        % edges2pixels
        selfEdgePixelSet(i) = edges2pixels(selfEdgeIDs(i),2); 
        edges2pixels(selfEdgeIDs(i),2) = 0;  % set the self edge 'pixel' to zero
    end
    % from edges2pixels, remove the rows who's second column has a zero
    edges2pixels = edges2pixels((edges2pixels(:,2)~=0),:);
    % nodeEdges may contain zeros for edgeIDs among nonzero entries. get rid of
    % the zeros
    numNodes = size(nodeEdges,1);
    for i=1:numNodes
       nodeEdgesList_i = nodeEdges(i,(nodeEdges(i,:)>0)); 
       numEdges = numel(nodeEdgesList_i);
       for j=1:numEdges
          nodeEdges2(i,j) = nodeEdgesList_i(j); 
       end
    end
    nodeEdges = nodeEdges2;
    clear nodeEdges2;
    % Now, after removing the self edges, the graph contains some junctions
    % with only two edges connecting to them. i.e. they are not junctions
    % anymore but just edges. At the moment we just keep them. and treat them
    % as 2 edge junctions.
end

if(saveMatrices)
    save('edges2pixels.mat','edges2pixels')
end


%% visualize graph
if(displayImg)
    [r,c] = ind2sub([sizeR sizeC],nodeInds);
    xy = [c r];
    % figure;gplot(adjacencyMat,xy,'-*');
    figure;gplotwl(adjacencyMat,xy);
    set(gca,'YDir','reverse');
    axis square
end
% n = size(adjacencyMat,1);
% k = 1:n;
% figure;gplot(adjacencyMat(k,k),xy(k,:),'-*')
% set(gca,'YDir','reverse');
% axis square

