function edgePolarities = getEdgePolarity(edgeListInds,edges2nodes,nodeListInd)

% returns edge polarities for each edge given in edgeListInds
% +1 for out
% -1 for in

% in edges2nodes, the 1st col is the start node, 2nd col is the end node.
% the edge is assumed to go from node1(out) to node2(in)

numEdges = numel(edgeListInds);

edgePolarities = zeros(numEdges,1);

for i=1:numEdges
    % for each edge, check the position of nodeListInd (1 or 2)
    nodes_i = edges2nodes(edgeListInds(i),:);
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