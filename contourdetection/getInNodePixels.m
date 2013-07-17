function inNodePixels = getInNodePixels(inEdgeIDs,nodeEdges,...
        nodeActivationVector,clusterNodeIDs)

% if all the edges connecting to a node are inEdges, the node is an inNode

% Inputs:
%   inEdges - edgeListInds for inactive edges inside active regions
%   nodeEdges - 
%   edgeIDsAll - 
%   nodeActivationVector - vector contains 1 for active nodes, 0 otherwise
%   clusterNodeIDs - 

% Output:
%   inNodePixels - nodeInds of nodes inside active regions.

inNodePixels = [];

% get inactive nodes
inactiveNodeListInds = find(nodeActivationVector<1);

% out of inactive nodes, get the ones inside active regions
for i=1:numel(inactiveNodeListInds)
    % get edges for this node
    nodeEdgeIDs_i = nodeEdges(inactiveNodeListInds(i),:);
    nodeEdgeIDs_i(1) = []; % first element is 
    nodeEdgeIDs_i = nodeEdgeIDs_i(nodeEdgeIDs_i>0);
    % if all the edges are inEdges, this node is an inNode
    numEdges_i = numel(nodeEdgeIDs_i);
    count = 0;      % count of inEdges at this node
    for j=1:numEdges_i
        x = find(inEdgeIDs==nodeEdgeIDs_i(j));
        if(~isempty(x))
            % this is an inEdge
            count = count + 1;
        end
    end
    if(count==numEdges_i)
        % this is an inNode
        % get the pixels
        pixInd1 = nodeEdges(inactiveNodeListInds(i),1);
        % check if it is a nested node
        r = find(clusterNodeIDs==pixInd1);
        if(~isempty(r))
            % is a nested node
            clusID = clusterNodeIDs(r,2);
            clusNodeListInds = find(clusterNodeIDs(:,2)==clusID);
            clusNodes = clusterNodeIDs(clusNodeListInds,1);
            inNodePixels = [inNodePixels; clusNodes];
        else
           inNodePixels = [inNodePixels; pixInd1]; 
        end
    end
    
    % get the pixels for this node, including clusterNodePixels
    
    % append it to the list of inNodePixels
end