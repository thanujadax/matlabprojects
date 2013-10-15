function [faceAdj,edges2cells,setOfCellsMat,listOfTwoRegionEdgeIDs,wsIDs] = getFaceAdjFromJnAdjGraph...
    (edgeIDs,nodeEdges,junctionTypeListInds,jAnglesAll_alpha,...
    boundaryEdgeIDs,edges2nodes,ws,edges2pixels)
% Inputs: adjacency graph of junctions (planar graph)
%   edgeIDs -
%   nodeEdges -
%   junctionTypeListInds -
%   jAnglesAll_alpha -
%   boundaryEdgeIDs - edges at the border of the image. Each of them belong
%   just to one cell
%   edges2nodes - 

% Outputs: 
%   faceAdj - adjacency graph of the faces of the input planar graph.
%   values are edgeIDs
%   edges2cells - each row corresponds to an edge. the two columns give you
%   the two cells that are connected by that edge. The order of edges
%   present is according to listOfEdgeIDs
%   setOfCellsMat - each row corresponds to a cell and contains the set of
%   edges bounding that cell as a row vector with zero padding.
%   listOfEdgeIDs - the edges considered in the faceAdjMatrix. Each of these edges
%   connects a pair of cells

MAX_EDGE_USAGE = 2;
MAX_BOUNDARY_EDGE_USAGE = 1;

% get edge list - keep a count of their usage i.e. the two regions each
% edge separates
numEdges = numel(edgeIDs);
edgeUsage = zeros(numEdges,2);  % col1: edgeID, col2: usage
edgeUsage(:,1) = edgeIDs;

% separately identify edges on the boundary. these edges will only bound
% one cell. p.s. there can be a few false negatives.

%setOfCells = []; % each row corresponds to a cell i.e. a set of edges enclosing a cell
cellInd = 0;
for i=1:numEdges
    
    % check usage
    currentEdgeID = edgeIDs(i);
    % start of debug code
    if(currentEdgeID==1051)
        aa = 00;
    end
    % end of debug code
        
    currentEdgeUsage = edgeUsage(i,2); 
    % if boundaryEdge, max usage is 1
    % check if the edge is a boundary edge.
    currentEdgeIsBoundary = 0;
    if(max(boundaryEdgeIDs==currentEdgeID))
        % is a boundary edge
        currentEdgeIsBoundary = 1;
        if(currentEdgeUsage>=MAX_BOUNDARY_EDGE_USAGE)
            continue
        end 
        % we don't want to initialize a loop with a boundary edge. - ????
%         continue
    else
        % not a boundary edge
        if(currentEdgeUsage>=MAX_EDGE_USAGE)
            continue
        end
    end
    % get next edge to complete clockwise loop at both ends
    currentNodeListInds = edges2nodes(i,:);    % row vector containing nodeListInds
    % check usage of next edge (both ends)
    nextEdgeIDs_2 = zeros(1,2);
    [nextEdgeIDs_2(1),~] = getNextEdge(currentEdgeID,currentNodeListInds(1),nodeEdges...
        ,junctionTypeListInds,jAnglesAll_alpha,edges2nodes,edgeIDs);
    [nextEdgeIDs_2(2),~] = getNextEdge(currentEdgeID,currentNodeListInds(2),nodeEdges...
        ,junctionTypeListInds,jAnglesAll_alpha,edges2nodes,edgeIDs);
    
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
        % if the current edge is also a boundary edge, pick edge2
        if(currentEdgeIsBoundary)
            edge1ok = 0;
        else
            if(nextEdgeUsage_2(1)<MAX_BOUNDARY_EDGE_USAGE)
                edge1ok = 1;
                % flag set. get loop            
            else
                edge1ok = 0;
            end
        end
    else
        % nextEdge(1) is not a boundary edge (most likely)
        if(nextEdgeUsage_2(1)<MAX_EDGE_USAGE)
            edge1ok = 1;
            % flag set. get loop
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
                % flag set. get loop
            else
                edge2ok = 0;
            end
        else
        % nextEdge(2) is not a boundary edge (most likely)
            if(nextEdgeUsage_2(2)<MAX_EDGE_USAGE)
                edge2ok = 1;
                % flag set. get loop.
            else
                edge2ok = 0;
            end
        end
    end   
    
    % find loops (closed contours)
    setOfEdges_loop = [];
    if(edge1ok)
        % look for loop containing edge1
        [setOfEdges_loop,edgeUsage_new] = getEdgeLoop(currentNodeListInds(1),...
            currentEdgeID,nodeEdges,junctionTypeListInds,jAnglesAll_alpha,edgeUsage,...
            boundaryEdgeIDs,MAX_EDGE_USAGE,MAX_BOUNDARY_EDGE_USAGE,...
            edges2nodes,edgeIDs);

    elseif(edge2ok)
        % look for loop containing edge2
        [setOfEdges_loop,edgeUsage_new] = getEdgeLoop(currentNodeListInds(2),...
            currentEdgeID,nodeEdges,junctionTypeListInds,jAnglesAll_alpha,edgeUsage,...
            boundaryEdgeIDs,MAX_EDGE_USAGE,MAX_BOUNDARY_EDGE_USAGE,...
            edges2nodes,edgeIDs);
    end
    
    if(~isempty(setOfEdges_loop) && setOfEdges_loop(1)~=0)
 
        cellInd = cellInd + 1;
        % setOfCells = [setOfCells; setOfEdges_loop];
        setOfCells{cellInd} = setOfEdges_loop;
        edgeUsage = edgeUsage_new;
    end
    
end
% create adjacency matrix for the cells. The coefficients correspond to the
% edgeID that connects the corresponding pair of cells
setOfCellsMat = setOfCells2Mat(setOfCells);

wsIDs = getWsIDsForCellIDs(ws,setOfCellsMat,edges2pixels,nodeEdges(:,1),...
            edges2nodes,edgeIDs);
        
% look for duplicate wsIDs
uniqueIDs = unique(wsIDs);
n = histc(wsIDs,uniqueIDs);
duplicateWsIDs = uniqueIDs(n>1);
listOfCellsToRemove = [];
[sizeR,sizeC] = size(ws);
while(~isempty(duplicateWsIDs))
    duplicateCellIDs_i = find(wsIDs==duplicateWsIDs(1));
    % find which cell is the biggest
    % get the internal pixels for each cell
    numDupCells_i = numel(duplicateCellIDs_i);
    numPixPerDupCell = zeros(numDupCells_i,1);
    for i=1:numDupCells_i
        boundaryPixels_ii = getBoundaryPixelsForCell(setOfCellsMat((duplicateCellIDs_i(i)),:)...
            ,edges2pixels,nodeEdges(:,1),edges2nodes,edgeIDs);
        [internalx,internaly] = getInternelPixelsFromBoundary(boundaryPixels_ii,sizeR,sizeC);
        % get the most popular label for the internal pixels
        internalPixels_ii = sub2ind([sizeR sizeC],internaly,internalx);
        numPixPerDupCell(i) = numel(internalPixels_ii);       
    end
    [~,pos_i] = max(numPixPerDupCell);
    cellToRem_i = duplicateCellIDs_i(pos_i);
    listOfCellsToRemove = [listOfCellsToRemove; cellToRem_i];
    duplicateWsIDs(1) = [];     % remove this wsID from the duplicates list
end

% removing detected duplicates
% sort cell ids in descending order
listOfCellsToRemove = sort(listOfCellsToRemove,'descend');
numCellsToRemove = numel(listOfCellsToRemove);
if(~isempty(listOfCellsToRemove))
    for i=1:numCellsToRemove
        setOfCellsMat(listOfCellsToRemove(i),:)=[];
        wsIDs(i) = [];
        
    end
end

numCells = size(setOfCellsMat,1);
cellList = 1:numCells; % row vector
% add the cellList (index) as the first col of setOfCells. This is done so
% that we can reuse getAdjacencyMat() to creage faceAdj.
setOfCellsMat_2 = [cellList' setOfCellsMat];
[faceAdj,edges2cells,~,listOfTwoRegionEdgeIDs] = getAdjacencyMat(setOfCellsMat_2);