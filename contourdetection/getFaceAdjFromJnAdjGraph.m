function faceAdj = getFaceAdjFromJnAdjGraph(A,edgeIDs,nodeEdges,junctionTypeListInds,...
    jAnglesAll_alpha,boundaryEdgeIDs,edges2nodes)
% input: adjacency graph of junctions (planar graph)
% output: adjacency graph of the faces of the input planar graph

MAX_EDGE_USAGE = 2;
MAX_BOUNDARY_EDGE_USAGE = 1;

% get edge list - keep a count of their usage i.e. the two regions each
% edge separates
numEdges = numel(edgeIDs);
edgeUsage = zeros(numEdges,2);  % col1: edgeID, col2: usage
edgeUsage(:,1) = edgeIDs;

% separately identify edges on the boundary. these edges will only bound
% one cell . p.s. there can be a few false negatives.

for i=1:numEdges
    % check usage
    currentEdgeID = edgeIDs(i);
    currentEdgeUsage = edgeUsage(i,2); 
    % if boundaryEdge, max usage is 1
    if(max(boundaryEdgeIDs==currentEdgeID))
        % is a boundary edge
        if(currentEdgeUsage>=MAX_BOUNDARY_EDGE_USAGE)
            continue
        end     
    else
        % most likely, not a boundary edge
        if(currentEdgeUsage>=MAX_EDGE_USAGE)
            continue
        end
    end
    % get next edge to complete clockwise loop at both ends
    currentNodeListInds = edges2nodes(i,:);    % row vector containing nodeListInds
    % check usage of next edge (both ends)
    nextEdgeIDs_2 = zeros(1,2);
    nextEdgeIDs_2(1) = getNextEdge(currentEdgeID,currentNodeListInds(1),nodeEdges...
        ,junctionTypeListInds,jAnglesAll_alpha);
    nextEdgeIDs_2(2) = getNextEdge(currentEdgeID,currentNodeListInds(2),nodeEdges...
        ,junctionTypeListInds,jAnglesAll_alpha);
    nextEdgeUsage_2 = zeros(1,2);
    nextEdgeUsage_2(1) = edgeUsage((edgeUsage(:,1)==nextEdgeIDs_2(1)),2);
    nextEdgeUsage_2(2) = edgeUsage((edgeUsage(:,1)==nextEdgeIDs_2(2)),2);
    % check if any of the candidates are ok in terms of usage
    edge1ok = 0;
    edge2ok = 0;
    % if at least one is ok
        % set usage for this edge and next edge
        % continue loop aggregation along next edge
    % else continue
    if(max(boundaryEdgeIDs==nextEdgeIDs_2(1)))
        % nextEdge(1) is a boundary edge
        if(nextEdgeUsage_2(1)<MAX_BOUNDARY_EDGE_USAGE)
            edge1ok = 1;
            % TODO: get loop
        else
            edge1ok = 0;
        end
    else
        % nextEdge(1) is not a boundary edge (most likely)
        if(nextEdgeUsage_2(1)<MAX_EDGE_USAGE)
            edge1ok = 1;
            % TODO: get loop
        else
            edge1ok = 0;
        end
    end
    if(~edge1ok)
        % if nextEdge(1) is not ok, check if nextEdge(2) is ok
        if(max(boundaryEdgeIDs==nextEdgeIDs_2(2)))
        % nextEdge(2) is a boundary edge
            if(nextEdgeUsage_2(2)<MAX_BOUNDARY_EDGE_USAGE)
                edge2ok = 1;
                % TODO: get loop
            else
                edge2ok = 0;
            end
        else
        % nextEdge(2) is not a boundary edge (most likely)
            if(nextEdgeUsage_2(2)<MAX_EDGE_USAGE)
                edge2ok = 1;
                % TODO: get loop
            else
                edge2ok = 0;
            end
        end
    end   
end
