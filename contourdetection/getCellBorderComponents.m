function [c_cellBorderEdgeIDs,c_cellBorderNodeIDs] = getCellBorderComponents(onEdgeInd,...
        edges2nodes,nodeEdges,edgeListInds)
    
% inputs:
%   edges2nodes - row index is edgeListInd
    
usedEdgeInds = [];

numEdges = size(onEdgeInd,1);
k = 0;
for i=1:numEdges
    % check if this edge is already accounted for
    edgeInd_i = onEdgeInd(i);
    if(sum(usedEdgeInds==edgeInd_i)==0)
        % edge_i is not yet used
        k = k+1;
        %usedEdgeInds(end+1) = edgeInd_i; % edge_i marked as used
        [c_cellBorderEdgeIDs{k},c_cellBorderNodeIDs{k}] = getBorderForCell(edgeInd_i,...
            onEdgeInd,edges2nodes,nodeEdges,edgeListInds);
        usedEdgeInds = [usedEdgeInds; c_cellBorderEdgeIDs{k}];
    end
end


function [visitedEdges, visitedNodes] = getBorderForCell(edgeInd_i,onEdgeInd,...
            edges2nodes,nodeEdges,edgeListInds)

% c_cellBorderPixels = [];
visitedNodes = [];  % records all the nodes visited already for this cell
visitedEdges = [];  % records all the edges visited already for this cell

nodesForEdge_i = edges2nodes(edgeInd_i,:);
node_i = nodesForEdge_i(1); % pick the first node to start traversing the cell
   
while(sum(visitedEdges==edgeInd_i)==0 && numel(node_i)==1)
    % add the pixels of edge_i to the cellBorderPixels list
    visitedEdges = [visitedEdges; edgeInd_i];
    
    % add the node pixel to the cellBorderPixels list
    visitedNodes(end+1) = node_i;
    
    % pick the next edge
    edgeIDs_for_node_i = nodeEdges(node_i,:);
    edgeIDs_for_node_i(1) = []; % first element is the nodePixInd

    [~,edgeInds_for_node_i] = ismember(edgeIDs_for_node_i,edgeListInds);
    active_edges_node_i = intersect(edgeInds_for_node_i,onEdgeInd);
    nextEdges = (active_edges_node_i==edgeInd_i);
    nextEdges = ~nextEdges;
    edgeInd_i = active_edges_node_i(nextEdges);
%     if(edgeInd_i==158)
%         a = 00;
%     end
    % pick the next node
    nextNodes = edges2nodes(edgeInd_i,:);
    nextNode_listInd = (nextNodes==node_i);
    nextNode_listInd = ~nextNode_listInd;
    node_i = nextNodes(nextNode_listInd);
    % start debug
    if(isempty(edgeInd_i))
        aa = 99;
    end
    % end debug
end
