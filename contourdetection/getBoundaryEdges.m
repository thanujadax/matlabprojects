function boundaryEdges = getBoundaryEdges(wsgraph,marginSize,ws_edgePixels,...
    edges2nodes,nodeEdges,edgeIDs)
% an additional thick margin is already added along the boundary of the
% image. Once the watershed oversegmentation is performed 
% Inputs:
%   wsgraph - graph obtained from oversegmenting edges from OFR using WS
%   barLength - 
%   ws_EdgePixels - each row contains the set of pixels for each edge (zero
%   padded at the end). no absolute edgeIDs.
%   edges2nodes - each row contains the nodeListInds of the two ends for
%   the edges in the order given in ws_EdgePixels
%   nodeEdges - each row contains the set of edgeIDs connected to each
%   node. the first col contains the pixelIDs of each node
%   edgeIDs - list of edgeIDs in the order given in ws_edgePixels

% Output
%   boundaryEdges - listInd (wrt the order in ws_edgePixels) of boundary
%   edges

visualize = 1;

[sizeR,sizeC] = size(wsgraph);
%% first pass
% get all possible boundary pixels from OFR
% separately for each boundary
topPixInd_col = find(wsgraph(marginSize,:)>0);
numTopPix = numel(topPixInd_col);
topPixInd_row = ones(numTopPix,1) .* marginSize;
topPixInd = sub2ind([sizeR sizeC],topPixInd_row,topPixInd_col');

botPixInd_col = find(wsgraph((sizeR-marginSize),:)>0);
numBotPix = numel(botPixInd_col);
botPixInd_row = ones(numBotPix,1) .* (sizeR-marginSize);
botPixInd = sub2ind([sizeR sizeC],botPixInd_row,botPixInd_col');

leftPixInd_row = find(wsgraph(:,marginSize)>0);
numLeftPix = numel(leftPixInd_row);
leftPixInd_col = ones(numLeftPix,1) .* marginSize;
leftPixInd = sub2ind([sizeR sizeC],leftPixInd_row,leftPixInd_col);

rightPixInd_row = find(wsgraph(:,(sizeC-marginSize)));
numRightPix = numel(rightPixInd_row);
rightPixInd_col = ones(numRightPix,1) .* (sizeC-marginSize);
rightPixInd = sub2ind([sizeR sizeC],rightPixInd_row,rightPixInd_col);

% out of the all possible boundary pixels from OFR, extract the ones that
% correspond to WS edges
edgeListInd_top = getEdgeListIndForPixInds(topPixInd,ws_edgePixels);
edgeListInd_bot = getEdgeListIndForPixInds(botPixInd,ws_edgePixels);
edgeListInd_left = getEdgeListIndForPixInds(leftPixInd,ws_edgePixels);
edgeListInd_right = getEdgeListIndForPixInds(rightPixInd,ws_edgePixels);

% first estimate of boundary edges (from above)
boundaryEdges_listInds = [edgeListInd_top'; edgeListInd_bot'; edgeListInd_left'; edgeListInd_right'];
boundaryEdges_listInds = unique(boundaryEdges_listInds); % contains the edgeListInds
%% 2nd pass
% if an edge is connected to boundary edges at both its ends, then this edge
% is also a boundary edge

% get the edges which are classified as 'not-boundary', so far
listAllEdges = 1:numel(edgeIDs);
nonClassifiedListInd = setdiff(listAllEdges,boundaryEdges_listInds);

% check each of them to see if anyone is connected to boundary edges at
% both ends. These are boundary edges.
newBoundaryEdges_listInds = [];
numNonClas = numel(nonClassifiedListInd);
for i = 1: numNonClas
    % get ends
    endsVec = edges2nodes(nonClassifiedListInd(i),:);
    % get all edges connected to each end and check if at least one of
    % those is a boundary edge
    edges_node1 = nodeEdges(endsVec(1),2:size(nodeEdges,2));
    edges_node1 = edges_node1(edges_node1>0);
    num_node1_edges = numel(edges_node1);
    for j=1:num_node1_edges
        edgeListInd_j = find(edgeIDs==edges_node1(j));
        isBoundaryEdge = find(boundaryEdges_listInds==edgeListInd_j);
        if(~isempty(isBoundaryEdge))
            % end1 has a boundary edge. check the other end
            edges_node2 = nodeEdges(endsVec(2),2:size(nodeEdges,2));
            edges_node2 = edges_node2(edges_node2>0);
            num_node2_edges = numel(edges_node2);
            for k=1:num_node2_edges
                edgeListInd_k = find(edgeIDs==edges_node2(k));
                clear isBoundaryEdge
                isBoundaryEdge = find(boundaryEdges_listInds==edgeListInd_k);
                if(~isempty(isBoundaryEdge))
                    % this edge is a boundary edge
                    newBoundaryEdges_listInds = ...
                        [newBoundaryEdges_listInds; nonClassifiedListInd(i)];
                end
            end
        end
    end
    % if both ends give 2 'yes's, classify this as a boundary edge
end
boundaryEdges_listInds = [boundaryEdges_listInds; newBoundaryEdges_listInds];
boundaryEdges_listInds = unique(boundaryEdges_listInds);
% TODO: get the actual edgeIDs for the boundary edges from the list inds
numBoundaryEdges = numel(boundaryEdges_listInds);
boundaryEdges = zeros(numBoundaryEdges,1);
for i=1:numBoundaryEdges
    boundaryEdges(i) = edgeIDs(boundaryEdges_listInds(i));
end

%% visualize
if(visualize)
    % visualize detected boundary edges
    imgTmp = zeros(sizeR,sizeC);
    % color all edge pixels : 1
    edgepix_all = ws_edgePixels(ws_edgePixels>0);
    imgTmp(edgepix_all) = 1;
    % color boundary edges with 0.5
    for i=1:numBoundaryEdges
        edgepix_i = ws_edgePixels(boundaryEdges_listInds(i),:);
        edgepix_i = edgepix_i(edgepix_i>0);
        imgTmp(edgepix_i) = 0.5;
    end
    figure;imagesc(imgTmp);title('boundary edges');
end
