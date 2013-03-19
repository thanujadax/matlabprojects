% ILP script 1

% max vote response image of the orientation filters
imFilePath = 'testMem4_V.png';
% votes for each orientation for each edge
load('orientedScoreSpace3D.mat') % loads the orientation filter scores
angleStep = 10; % 10 degrees discretization step of orientations

% param
cNode = 1;          % scaling factor for the node cost coming from gaussian normal distr.
sig = 45;          % standard deviation(degrees) for the node cost function's gaussian distr.
midPoint = 180;     % angle difference of an edge pair (in degrees) for maximum cost 

imIn = imread(imFilePath);
% watershed segmentation
ws = watershed(imIn);
[sizeR,sizeC] = size(ws);
%% generate graph from the watershed edges
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels] = getGraphFromWS(ws);

%% Edge priors
% edge priors - from orientation filters
edgePriors = getEdgePriors(orientedScoreSpace3D,edges2pixels);

%% Edge pairs - Junction costs
nodeInds = nodeEdges(:,1);
j4ListInd = find((nodeEdges(:,5))>0);       % indices of J4 in order
j3ListInd = find((nodeEdges(:,5))==0);      % indices of J3 in order
numJ4 = numel(j4ListInd);
numJ3 = numel(j3ListInd);
% J3
% j3Edges = zeros(numJ3,3);
j3Edges = nodeEdges(j3ListInd,2:4);             % order of indices given by j3Ind
% angles
j3Angles = getNodeAngles(j3ListInd,nodeInds,j3Edges,edges2pixels,orientedScoreSpace3D,...
                            sizeR,sizeC,angleStep);                        
j3dTheta = getAngleDifferences(j3Angles);
% calculate the cost for the angle differences
j3NodeAngleCost = getNodeAngleCost(j3dTheta,midPoint,sig,cNode);

% J4
% j4Edges = zeros(numJ4,4);
j4Edges = nodeEdges(j4ListInd,2:5);             % order of indices given by j4Ind
% angles
j4Angles = getNodeAngles(j4ListInd,nodeInds,j4Edges,edges2pixels,orientedScoreSpace3D,...
                            sizeR,sizeC,angleStep);
j4dTheta = getAngleDifferences(j4Angles);
% calculate the cost for the angle differences
j4NodeAngleCost = getNodeAngleCost(j4dTheta,midPoint,sig,cNode);

%% ILP
% cost function to minimize
% state vector x: {edges*2}{J3*4}{J4*7}
numEdges = size(edges2nodes,1);
numJunctions = numJ3 + numJ4;
% tot num of int variables = 2*numEdges + 4*numJ3 + 7*numJ4
% coeff (unary prior) for turning off each edge = +edgePriors (col vector)
% coeff (unary prior) for turning on each edge = -edgePriors (col vector)
% coeff for turning off J3s: min(j3NodeAngleCost): max(j3NodeAngleCost)
% coeff for turning on J3-config(1 to 3): j3NodeAngleCost
% coeff for turning off J4s: max(j3NodeAngleCost)
% coeff for turning on J4-config(1 to 6): j4NodeAngleCost
f = getILPcoefficientVector(edgePriors,j3NodeAngleCost,j4NodeAngleCost);
% constraints
% equality constraints and closedness constrains in Aeq matrix
[Aeq,beq] = getEqConstraints(numEdges,j3Edges,j4Edges);
% ILP
x0 = getInitValues(numEdges,numJ3,numJ4);    