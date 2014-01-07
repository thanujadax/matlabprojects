function [nodeEdges,nodeIndsNoDuplicates] = getNodeEdges(nodeInds,edgePixLabels,connectedJunctionIDs,sizeR,sizeC)
% Inputs:
%   nodeInd - array of junction indices
%   edgePixLabels - N-by-2 array of edge labels for each pixel, given by
%   the index wrt the original image (watershed)
%   connectedJunctionIDs - list of clusterd junction indices with the associated
%   cluster label.

% Output:
%   nodeEdges - array with edge labels corresponding to each junction. each
%   row -> jn,edge1,edge2,edge3,edge4
%       edgeIDs sorted in ascending order 20140107
%   nodeIndsNoDuplicates - list of the pixel indices of the nodes provided
%   in nodeEdges. For the detected clustered nodes in connectedJunctions 

% for each node, get the neighbors
% get the edgeID of the neighbors
if(size(connectedJunctionIDs,2)==2)
    numClusters = max(connectedJunctionIDs(:,2));
    numClusteredNodes = size(connectedJunctionIDs,1);
else
    numClusters = 0;
    numClusteredNodes = 0;
end
% number of nodes after combining clusters
numNodes0 = numel(nodeInds);
numNodesCombined = numNodes0 - numClusteredNodes + numClusters;

% create a list of nodeInds without duplicates
nodeIndsNoDuplicates = nodeInds;
for i=1:numClusters
    cNodesListInd = find(connectedJunctionIDs(:,2)==i);
    cNodes_i = connectedJunctionIDs(cNodesListInd,1);
    numCnodes_i = numel(cNodes_i);
    for j=2:numCnodes_i
        nodeIndsNoDuplicates = nodeIndsNoDuplicates...
                    (nodeIndsNoDuplicates(:,1)~=cNodes_i(j),:);
    end    
end

nodeEdges = zeros(numNodesCombined,5);

for i=1:numNodesCombined
    thisNodeIndex = nodeIndsNoDuplicates(i);
    neighborInd = getNeighbors(thisNodeIndex,sizeR,sizeC);
    numNeighbors = numel(neighborInd);
    nodeEdges(i,1) = thisNodeIndex;
    k = 1;
    for j=1:numNeighbors
        neighborListInd = find(edgePixLabels(:,1)==neighborInd(j));
        if(~isempty(neighborListInd))
            % there's an edge for this neighbor
            k = k + 1;
            edgeId = edgePixLabels(neighborListInd,2);
            nodeEdges(i,k) = edgeId;
        end
    end
    % check if it has directly neighboring nodes
    connectedJInd = find(connectedJunctionIDs(:,1)==thisNodeIndex);
    if(~isempty(connectedJInd))
        % this junction has other directly neighboring junctions
        duplicateLabel = connectedJunctionIDs(connectedJInd,2);
        neighboringCJListInd = find(connectedJunctionIDs(:,2)==duplicateLabel(1));
        
        numNeighJ = numel(neighboringCJListInd);
        for m=2:numNeighJ      % start with 2 to skip the current junction node
            % for each neighboring junction,
            neighJ_ind = connectedJunctionIDs(neighboringCJListInd(m),1);
            % get its edges and,
            % add its edges to the list of edges under 'thisNodeIndex'
            jNeighborInd = getNeighbors(neighJ_ind,sizeR,sizeC);
            jNeighborInd = jNeighborInd(jNeighborInd~=thisNodeIndex);
            numJNeighbors = numel(jNeighborInd);
            for j=1:numJNeighbors
                jNeighborListInd = find(edgePixLabels(:,1)==jNeighborInd(j));
                if(~isempty(jNeighborListInd))
                    % there's an edge for this neighbor
                    % check first, if this edge is already in the list for
                    % thisNodeIndex
                    edgeId = edgePixLabels(jNeighborListInd,2);
                    if(isempty(find(nodeEdges(i,:)==edgeId)) )
                        k = k + 1;
                        nodeEdges(i,k) = edgeId;
                    end
                end
            end                        
        end
    end
    % sort edgeIDs in ascending order
    if(k>1)
        edgeList_i = nodeEdges(i,2:k);
        edgeList_i = sort(edgeList_i);
        nodeEdges(i,2:k) = edgeList_i;
    end
end