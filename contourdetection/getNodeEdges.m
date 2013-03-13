function nodeEdges = getNodeEdges(nodeInds,edgePixLabels,connectedJunctionIDs,sizeR,sizeC)
% Inputs:
%   nodeInd - array of junction indices
%   edgePixLabels - N-by-2 array of edge labels for each pixel, given by
%   the index wrt the original image (watershed)

% Output:
%   nodeEdges - array with edge labels corresponding to each junction. each
%   row -> jn,edge1,edge2,edge3,edge4

% for each node, get the neighbors
% get the edgeID of the neighbors
numClusters = max(connectedJunctionIDs(:,2));
numClusteredNodes = size(connectedJunctionIDs,1);
% number of nodes after combining clusters
numNodes0 = numel(nodeInds);
numNodesCombined = numNodes0 - numClusteredNodes + numClusters;

nodeEdges = zeros(numNodesCombined,5);
% nodeEdges(:,1) = nodeInds;

for i=1:numel(nodeInds)
    thisNodeIndex = nodeInds(i);
    neighborInd = getNeighbors(thisNodeIndex,sizeR,sizeC);
    numNeighbors = numel(neighborInd);
    nodeEdges(i,1) = thisNodeIndex;
    k = 1;
    for j=1:numNeighbors
        neighborListInd = find(edgePixLabels(:,1)==neighborInd(j));
        if(~isempty(neighborListInd))
            % there's an edge for this neighbor
            k = k + 1;
            nodeEdges(i,k) = edgePixLabels(neighborListInd,2);
        end
    end
    % check if it has directly neighboring nodes
    connectedJInd = find(connectedJunctionIDs(:,1)==thisNodeIndex);
    if(~isempty(connectedJInd))
        % this junction has other directly neighboring junctions
        duplicateLabel = connectedJunctionIDs(connectedJInd,2);
        neighboringCJListInd = find(connectedJunctionIDs(:,2)==duplicateLabel);
        
        numNeighJ = numel(neighboringCJListInd);
        for m=1:numNeighJ
            % for each neighboring junction,
            neighJ_ind = connectedJunctionIDs(neighboringCJListInd(m),1);
            % get its edges and,
            % add its edges to the list of edges under 'thisNodeIndex'
            jNeighborInd = getNeighbors(neighJ_ind,sizeR,sizeC);
            numJNeighbors = numel(jNeighborInd);
            for j=1:numJNeighbors
                jNeighborListInd = find(edgePixLabels(:,1)==jNeighborInd(j));
                if(~isempty(jNeighborListInd))
                    % there's an edge for this neighbor
                    k = k + 1;
                    nodeEdges(i,k) = edgePixLabels(jNeighborListInd,2);
                end
                % remove the duplicate junction from nodeInds
                % b = a(a~=3);
                nodeInds = nodeInds(nodeInds(:,1)~=jNeighborInd(j),:);
            end
            
            
        end
    end
end