function edgePolarities = getEdgePolarity(edgeListInds,edges2nodes,nodeListInd,...
                nodeEdgeLIDs)

% returns edge polarities for each edge given in edgeListInds, wrt to the
% assigned direction for the edges according to edges2nodes.

% +1 for out
% -1 for in

% in edges2nodes, the 1st col is the start node, 2nd col is the end node.
% the edge is assumed to go from node1(out) to node2(in)

numNodeEdges = numel(nodeEdgeLIDs);

edgePolarities = zeros(numNodeEdges,1);

for i=1:numNodeEdges
    % for each edge, check the position of nodeListInd (1 or 2)
    nodes_i = edges2nodes(nodeEdgeLIDs(i),:);
    if(nodes_i(1)==nodeListInd)
        % outwards
        edgePolarities(i) = 1;
    elseif(nodes_i(2)==nodeListInd)
        % inwards
        edgePolarities(i) = -1;        
    else
        disp('ERROR: getEdgePolarity.m: matching node not found')
    end
end