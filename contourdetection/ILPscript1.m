% ILP script 1

% max vote response image of the orientation filters
imFilePath = 'testMem4_V.png';
% votes for each orientation for each edge
load('orientedScoreSpace3D.mat') % loads the orientation filter scores

imIn = imread(imFilePath);
% watershed segmentation
ws = watershed(imIn);
% generate graph from the watershed edges
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels] = getGraphFromWS(ws);

% state vector x: {edges}{J3}{J4}
numEdges = size(edges2nodes,1);
numJunctions = size(nodeEdges,1);

% priors
% edge priors - from orientation filters
edgePriors = getEdgePriors(orientedScoreSpace3D,edges2pixels);
% J
nodeInds = nodeEdges(:,1);
j4Ind = nodeInds((nodeEdges(:,5))>0);
j3Ind = nodeInds((nodeEdges(:,5))==0);
numJ4 = numel(j4Ind);
numJ3 = numel(j3Ind);
% J3
j3edges = zeros(numJ3,3);
j3edges = nodeEdges(j3Ind,2:4);
% angles

% J4
j4edges = zeros(numJ4,4);
j4edges = nodeEdges(j4Ind,2:5);
% angles

