function faceAdj = getFaceAdjFromJnAdjGraph(A,edgeIDs,nodeEdges)
% input: adjacency graph of junctions (planar graph)
% output: adjacency graph of the faces of the input planar graph

% draw input graph in XY plane
% for each node, get (x,y) coordinates

% traverse the graph to pick one cell at a time

% get edge list - keep a count of their usage i.e. the two regions each
% edge separates
numEdges = numel(edgeIDs);

% create region list - each region stores its bounding set of edges

numIter = numEdges * 2;