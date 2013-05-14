function jNodeListInds = getJunctionsForEdges(edges2nodes,edgeListInds)

% edges2nodes: contains the 2 nodes corresponding to the edge. The list of
% edges is according the usage in the cost function (doesn't have to be the
% same as the set of edgeIDs)
% nodeListInd: list of junction node indices in the order they are used in
% the cost function
% edgeListInds: the list indices (according to the cost function) of the
% edges whose corresponding junction nodes have to be returned

numEdgesIn = numel(edgeListInds);
clear jNodeListInds;
jNodeListInds = unique(edges2nodes(edgeListInds,:));

