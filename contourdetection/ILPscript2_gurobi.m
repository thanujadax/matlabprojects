% ILP script 2

isToyProb = 0;
useGurobi = 1;

% max vote response image of the orientation filters
if(isToyProb)
    imFilePath = 'testMem4_V.png';
    % votes for each orientation for each edge
    load('orientedScoreSpace3D.mat') % loads the orientation filter scores
else
    imFilePath = 'stem_256x_t02_V.png';
    % votes for each orientation for each edge
    load('orientedScoreSpace3D_stem256x.mat') % loads the orientation filter scores
end

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
disp('creating graph from watershed boundaries...');
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs] = getGraphFromWS(ws);
nodeInds = nodeEdges(:,1);                  % indices of the junction nodes
j4ListInd = find((nodeEdges(:,5))>0);       % nodeInds list indices of J4 in order
j3ListInd = find((nodeEdges(:,5))==0);      % nodeInds list indices of J3 in order
clusterNodeIDs = connectedJunctionIDs(:,1); % indices of the clustered junction nodes
disp('graph created!')
wsBoundariesFromGraph = zeros(sizeR,sizeC);
wsBoundariesFromGraph(nodeInds(j3ListInd)) = 0.7;    % j3 nodes
wsBoundariesFromGraph(nodeInds(j4ListInd)) = 0.4;    % j4 nodes
wsBoundariesFromGraph(clusterNodeIDs) = 0.5;    % j4 nodes
wsBoundariesFromGraph(edges2pixels(edges2pixels>0)) = 1; % edge pixels
figure;imagesc(wsBoundariesFromGraph);title('boundaries from graph') 
disp('preparing coefficients for ILP solver...')
%% Edge priors
% edge priors - from orientation filters
edgePriors = getEdgePriors(orientedScoreSpace3D,edges2pixels);

%% Edge pairs - Junction costs
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
%% solver
if(useGurobi)
    disp('using Gurobi ILP solver...');
    model.A = sparse(Aeq);
    model.rhs = beq;
    model.obj = f';
    model.sense = '=';  % for the constraints given in A
    model.vtype = 'B';  % binary variables
    model.modelname = 'contourDetectionILP1';
    
    params.LogFile = 'gurobi.log';
    
    resultGurobi = gurobi(model,params);
    x = resultGurobi.x;
    
    
else
    % Matlab ILP solver
    disp('using MATLAB ILP solver...');
    % Initial values for the state variables
    % x0 = getInitValues(numEdges,numJ3,numJ4);  % TODO: infeasible!!
    numStates = size(f,1);
    maxIterationsILP = numStates * 1000000;
    options = optimoptions('bintprog','MaxIter',maxIterationsILP,...
                    'MaxTime',5000000);
    % options = struct('MaxTime', 5000000);
    disp('running ILP...');
    t1 = cputime;
    [x,fval,exitflag,optOutput] = bintprog(f,[],[],Aeq,beq,[],options);
    t2 = cputime;
    timetaken = t2-t1
end


%% visualize
% get active edges and active nodes from x
ilpSegmentation = zeros(sizeR,sizeC);
% active edges
onStateEdgeXind = 2:2:(numEdges*2);
onEdgeStates = x(onStateEdgeXind);
onEdgeInd = find(onEdgeStates==1);
onEdgePixelInds = getPixSetFromEdgeIDset(onEdgeInd,edges2pixels);
ilpSegmentation(onEdgePixelInds) = 1;

% active J3 nodes
offStateJ3Xind = (numEdges*2+1):4:(numEdges*2+numJ3*4-1);
offJ3PixStates = x(offStateJ3Xind);
onJ3PixInd = nodeInds(j3ListInd(offJ3PixStates==0));
ilpSegmentation(onJ3PixInd) = 0.7;
% active J4 nodes
offStateJ4Xind = (numEdges*2+numJ3*4+1):7:(numEdges*2+numJ3*4+numJ4*7-1);
offJ4PixStates = x(offStateJ4Xind);
onJ4PixInd = nodeInds(j4ListInd(offJ4PixStates==0));
ilpSegmentation(onJ4PixInd) = 0.3;
% active clustered nodes
activeNodesJ3J4 = [onJ3PixInd; onJ4PixInd];
numJ3J4Active = numel(activeNodesJ3J4);
activeClustNodeInd = 0;
for i=1:numJ3J4Active
    clustLabel = connectedJunctionIDs((connectedJunctionIDs(:,1)==activeNodesJ3J4(i)),2);
    clustNodeInd = connectedJunctionIDs((connectedJunctionIDs(:,2)==clustLabel),1);
    foundClustNode = find(activeNodesJ3J4==clustNodeInd);
    if(isempty(foundClustNode))
        activeClustNodeInd = [activeClustNodeInd; clustNodeInd];
    end
end
activeClustNodeInd = activeClustNodeInd(activeClustNodeInd~=0);
ilpSegmentation(activeClustNodeInd) = 0.5;
%onPixelInds = [onEdgePixelInds; onJ3PixInd; onJ4PixInd];
%ilpSegmentation(onPixelInds) = 1;
figure;imagesc(ilpSegmentation);title('ILP contours');
