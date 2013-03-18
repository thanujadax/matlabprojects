% ILP script 1

% max vote response image of the orientation filters
imFilePath = 'testMem4_V.png';
% votes for each orientation for each edge
load('orientedScoreSpace3D.mat') % loads the orientation filter scores
angleStep = 10; % 10 degrees discretization step of orientations

imIn = imread(imFilePath);
% watershed segmentation
ws = watershed(imIn);
[sizeR,sizeC] = size(ws);
%% generate graph from the watershed edges
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels] = getGraphFromWS(ws);

% state vector x: {edges}{J3}{J4}
numEdges = size(edges2nodes,1);
numJunctions = size(nodeEdges,1);

%% Edge priors
% edge priors - from orientation filters
edgePriors = getEdgePriors(orientedScoreSpace3D,edges2pixels);

%% Edge pairs - Junction costs
nodeInds = nodeEdges(:,1);
j4Ind = nodeInds((nodeEdges(:,5))>0);       % indices of J4 in order
j3Ind = nodeInds((nodeEdges(:,5))==0);      % indices of J3 in order
numJ4 = numel(j4Ind);
numJ3 = numel(j3Ind);
% J3
j3Edges = zeros(numJ3,3);
j3Edges = nodeEdges(j3Ind,2:4);             % order of indices given by j3Ind
% angles
j3Angles = getNodeAngles(j3Ind,j3Edges,edges2pixels,orientedScoreSpace3D,...
                            sizeR,sizeC,angleStep);

% J4
j4Edges = zeros(numJ4,4);
j4Edges = nodeEdges(j4Ind,2:5);             % order of indices given by j4Ind
% angles
j4Angles = getNodeAngles(j4Ind,j4Edges,edges2pixels,orientedScoreSpace3D,...
                            sizeR,sizeC,angleStep);



    