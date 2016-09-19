function [connectedJunctionIDs,pEdges2nodes] = getClusterNodesAndPsEdges(wsJ)

% Inputs:
%   wsJ - matrix containing 1 for junction pixels detected. others 0.

% Output:
%   connectedJunctionIDs - N-by-2 array. each cluster of junctions assigned
%   a unique label from 1 to n. n = number of clusters.
%   pEdges2nodes - psuedoEdges. contains the nodeInds which are being
%   connected


% extract edges with zero pixel length
% junction nodes (pixels) which are next to each other
% wsJ contains 1 for junction nodes, others zero.

[sizeR,sizeC] = size(wsJ);
ind4J = find(wsJ);

fourNH_J = zeros(sizeR,sizeC);
eightNH_J = zeros(sizeR,sizeC);
numJunctionPixels = numel(ind4J);  % watershed junction pixels (3 and 4 Js)
for i=1:numJunctionPixels
    % calculate n.o. 4 neighbors
    ind = ind4J(i);
    [r,c] = ind2sub([sizeR sizeC],ind);
    nh4 = zeros(1,4);
    nh8 = zeros(1,4);
    % N1
    if((r-1)>0)
        nh4(1) = wsJ(r-1,c);
    end
    % N2
    if((r+1)<=sizeR)
        nh4(2) = wsJ(r+1,c);
    end
    % N3
    if((c-1)>0)
        nh4(3) = wsJ(r,c-1);
    end
    % N4
    if((c+1)<=sizeC)
        nh4(4) = wsJ(r,c+1);
    end
    % N5
    if((r-1)>0 && (c-1)>0)
        nh8(1) = wsJ(r-1,c-1);
    end
    % N6
    if((r-1)>0 && (c+1)<=sizeC)
        nh8(2) = wsJ(r-1,c+1);
    end
    % N7
    if((r+1)<=sizeR && (c-1)>0)
        nh8(3) = wsJ(r+1,c-1);
    end
    % N8
    if((r+1)<=sizeR && (c+1)<=sizeC)
        nh8(4) = wsJ(r+1,c+1);
    end
    fourNH_J(ind) = sum(nh4);
    eightNH_J(ind) = sum(nh8);
end

indJClusterPixels = find(fourNH_J>0);  % gets junction pixels which has neighbors
% make a N-by-2 array of such junction pixels having immediate neighboring
% junctions
numJunctionClustPix = numel(indJClusterPixels);
if(numJunctionClustPix==0)
    % no clustered junction nodes
    connectedJunctionIDs = 0;
else
    connectedJunctionIDs = zeros(numJunctionClustPix,2);
    connectedJunctionIDs(:,1) = indJClusterPixels;
    junctionLabel = 0;
    for i=1:numJunctionClustPix
        % look for the neighbors and give them the same label
        jLabelCurrent = connectedJunctionIDs(i,2);
        if(jLabelCurrent==0)
            junctionLabel = junctionLabel + 1;
            % assign label to this junction and to its neighbors and its
            % neighbors neighbors
            junctionInd = connectedJunctionIDs(i,1);
            connectedJunctionIDs = labelJunctionCluster4Neighbors(junctionInd,connectedJunctionIDs,junctionLabel,...
                        sizeR,sizeC,fourNH_J);
        end
    end
    % connectedJunctionIDs contain the same ID for each node that is connected
    % together with zero length edges
    % nodeZeroEdges - store node - edge1,edge2 etc for these zero length edges
end

% pseudo edges between only-8-neighboring pixels
psuedoEdgeNodeInds = find(eightNH_J>0);
numPsuedoEdgeNodes = numel(psuedoEdgeNodeInds);
if(numPsuedoEdgeNodes==0)
    pEdges2nodes = 0;
else
    pEdges2nodes = zeros(1,2);
    k = 0;
    for i = 1: numPsuedoEdgeNodes
        % get 8 neighbors
        neighbors = get8NeighborsWithout4(psuedoEdgeNodeInds(i),sizeR,sizeC);
        neighbors = intersect(neighbors,psuedoEdgeNodeInds);
        if(numel(neighbors)>0)
            for j=1:numel(neighbors)
                % check if pair of nodes is in same cluster
                % check if psuedoEdge exists between them
                % if not assign pseudoedge between them
                if(isNotInCluster...
                        (psuedoEdgeNodeInds(i),neighbors(k),connectedJunctionIDs) ...
                        && noPsEdge(psuedoEdgeNodeInds(i),neighbors(k),pEdges2nodes))
                    k = k+1;
                    pEdges2nodes(k,1) = psuedoEdgeNodeInds(i);
                    pEdges2nodes(k,2) = neighbors(k);
                end
            end
        end
    end
end

function noEdge = noPsEdge(node1,node2,pEdges2nodes)

% check col1
node1ListInds_col1 = find(pEdges2nodes(:,1)==node1);
node1ListInds_col2 = find(pEdges2nodes(:,2)==node1);

noEdge = 1;

if(~isempty(node1ListInds_col1))
    if(sum(pEdges2nodes(node1ListInds_col1,2) == node2) > 0)
        noEdge = 0;
    end
end

if(~isempty(node1ListInds_col2))
    if(sum(pEdges2nodes(node1ListInds_col2,1) == node2) > 0)
        noEdge = 0;
    end
end

function notInCluster = isNotInCluster...
    (node1,node2,connectedJunctionIDs)

node1_clusListInd = find(connectedJunctionIDs(:,1)==node1);

if(~isempty(node1_clusListInd))
    node2_clusListInd = find(connectedJunctionIDs(:,1)==node2);
    if(~isempty(node2_clusListInd))
        node1_clusID = connectedJunctionIDs(node1_clusListInd,2);
        node2_clusID = connectedJunctionIDs(node2_clusListInd,2);
        if(node1_clusID==node2_clusID)
            notInCluster = 0;
        else
            notInCluster = 1;
        end
    else
        notInCluster = 1;
    end
else
    notInCluster = 1;
end
