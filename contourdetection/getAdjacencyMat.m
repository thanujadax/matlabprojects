function [adjacencyMat, edges2nodes,selfEdgeIDs,listOfEdgeIDs] = getAdjacencyMat(nodeEdges)
% Input:
%   nodeEdges: gives the list of edgeIDs connected to each junction node
%       each row is -> junctionInd,edge1, edge2, edge3, edge4, ..
[numNodes, numEdgesPerNode] = size(nodeEdges);
adjacencyMat = zeros(numNodes);

numEdges = max(max(nodeEdges(:,2:numEdgesPerNode)));
edges2nodes = zeros(numEdges,2);

sid = 0;
k = 0;
for i=1:numEdges
    % for each edge, find the two corresponding nodes at its ends
    [R,C] = find(nodeEdges(:,2:numEdgesPerNode)==i);
    % R has the list indices of the junctions corresponding to edge i
    if(numel(R)==2)
        % assign to adjacencyMat
        nodeInd = nodeEdges(R,1);
        j1 = find(nodeEdges(:,1)==nodeInd(1));
        j2 = find(nodeEdges(:,1)==nodeInd(2));        
        if(j1~=j2)
            % assign edgeId to the adjMat
            adjacencyMat(j1,j2) = i; 
            adjacencyMat(j2,j1) = i;
            % also add the entries to edges2nodes
            edges2nodes(i,1) = j1;
            edges2nodes(i,2) = j2;
            % also, add the edgeID to the listOfEdgeIDs
            k = k + 1;
            listOfEdgeIDs(k) = i;
        else
            sid = sid + 1;
            selfEdgeIDs(sid) = i;
        end
    elseif(numel(R)==1)
        % if 1, it contains a self edge.
        % disp('warning:getAdjacencyMat - edge skipped')       
        sid = sid + 1;
        selfEdgeIDs(sid) = i;
        
    else
        % disp('warning:getAdjacencyMat - edge skipped')
        % i
        sid = sid + 1;
        selfEdgeIDs(sid) = i;
        
    end
end

if(sid==0)
    selfEdgeIDs = 0;
end