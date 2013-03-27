% ILP script 3

isToyProb = 0;
useGurobi = 1;


orientations = 0:10:350;
barLength = 11; % should be odd
barWidth = 3; %
threshFrac = 0;
medianFilterH = 0;
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
cEdge = 10;
cNode = 100;          % scaling factor for the node cost coming from gaussian normal distr.
sig = 50;          % standard deviation(degrees) for the node cost function's gaussian distr.
midPoint = 180;     % angle difference of an edge pair (in degrees) for maximum cost 
% param for exp cost function
decayRate = 0.02;


imIn = imread(imFilePath);
% watershed segmentation
ws = watershed(imIn);
[sizeR,sizeC] = size(ws);
%% generate graph from the watershed edges
disp('creating graph from watershed boundaries...');
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs] = getGraphFromWS(ws);
nodeInds = nodeEdges(:,1);                  % indices of the junction nodes
junctionTypeListInds = getJunctionTypeListInds(nodeEdges);
% col1 has the listInds of type J2, col2: J3 etc. listInds are wrt
% nodeInds list of pixel indices of the detected junctions
clusterNodeIDs = connectedJunctionIDs(:,1); % indices of the clustered junction nodes
disp('graph created!')
wsBoundariesFromGraph = zeros(sizeR,sizeC);
wsBoundariesFromGraph(nodeInds) = 0.7;          % junction nodes
wsBoundariesFromGraph(clusterNodeIDs) = 0.5;    % cluster nodes
[nre,nce] = size(edges2pixels);  % first column is the edgeID
edgepixels = edges2pixels(:,2:nce);
wsBoundariesFromGraph(edgepixels(edgepixels>0)) = 1; % edge pixels
figure;imagesc(wsBoundariesFromGraph);title('boundaries from graph') 
disp('preparing coefficients for ILP solver...')
%% Edge priors
% edge priors - from orientation filters
edgePriors = getEdgePriors(orientedScoreSpace3D,edges2pixels);

%% Edge pairs - Junction costs
[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);
for i=1:numJtypes
    numJ(i) = sum(junctionTypeListInds(:,i)>0);
end

jEdges = getEdgesForAllNodeTypes(nodeEdges,junctionTypeListInds);
% jEdges{i} - cell array. each cell corresponds to the set of edges for the
% junction of type i (type1 = J2). A row of a cell corresponds to a node of
% that type of junction.
jAnglesAll = getNodeAnglesForAllJtypes(junctionTypeListInds,...
    nodeInds,jEdges,edges2pixels,orientedScoreSpace3D,sizeR,sizeC,angleStep);
% jAnglesAll{i} - cell array. each row of a cell corresponds to the set of angles for each
% edge at each junction of type 1 (= J2)

% angle differences for all edge combinations of all the junction types
dTheta = cell(1,numJtypes);
for i=1:numJtypes
    jAngles_i = jAnglesAll{i};
    dTheta{i} = getAngleDifferences(jAngles_i);
end
% angle costs
nodeAngleCosts = cell(1,numJtypes);
for i=1:numJtypes
    dTheta_i = dTheta{i};
    if(dTheta_i<0)
        % no such angles for this type of junction
    else
        %nodeAngleCosts{i} = getNodeAngleCost(dTheta_i,midPoint,sig,cNode);
        nodeAngleCosts{i} = getNodeAngleCostExp(dTheta_i,midPoint,decayRate,cNode);
    end
end


%% ILP
% cost function to minimize
% state vector x: {edges*2}{J3*4}{J4*7}
numEdges = size(edges2nodes,1);
numJunctions = numel(nodeInds);
% tot num of int variables = 2*numEdges + 4*numJ3 + 7*numJ4
% coeff (unary prior) for turning off each edge = +edgePriors (col vector)
% coeff (unary prior) for turning on each edge = -edgePriors (col vector)
% coeff for turning off J3s: min(j3NodeAngleCost): max(j3NodeAngleCost)
% coeff for turning on J3-config(1 to 3): j3NodeAngleCost
% coeff for turning off J4s: max(j3NodeAngleCost)
% coeff for turning on J4-config(1 to 6): j4NodeAngleCost

% f = getILPcoefficientVector(edgePriors,j3NodeAngleCost,j4NodeAngleCost);
scaledEdgePriors = edgePriors.*cEdge;
f = getILPcoefficientVector2(scaledEdgePriors,nodeAngleCosts);
% constraints
% equality constraints and closedness constrains in Aeq matrix
[Aeq,beq] = getEqConstraints2(numEdges,jEdges,edges2pixels);
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
    Initial values for the state variables
    x0 = getInitValues(numEdges,numJ3,numJ4);  % TODO: infeasible!!
    numStates = size(f,1);
    maxIterationsILP = numStates * 1000000;
    options = optimset('MaxIter',maxIterationsILP,...
                    'MaxTime',5000000);
    options = struct('MaxTime', 5000000);
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
% consider the edgeID given in the first col of edges2pixels?? no need for
% this since we use edgepixels array which is already sans the skipped
% edges
onStateEdgeXind = 2:2:(numEdges*2);
onEdgeStates = x(onStateEdgeXind);
onEdgeInd = find(onEdgeStates==1);
onEdgePixelInds = getPixSetFromEdgeIDset(onEdgeInd,edgepixels);
ilpSegmentation(onEdgePixelInds) = 1;
% active nodes
fIndStop = 2*numEdges;
for i=1:numJtypes
    % for each junction type
    % get the list of junctions and check their states in vector 'x'
    junctionNodesListListInds_i = find(junctionTypeListInds(:,i));
    if(~isempty(junctionNodesListListInds_i))
        junctionNodesListInds_i = junctionTypeListInds(junctionNodesListListInds_i,i);
        numJnodes_i = numel(junctionNodesListInds_i);
        % get the indices (wrt x) for the inactivation of the junctions
        numEdgePJ_i = i+1;
        numStatePJ_i = nchoosek(numEdgePJ_i,2)+1; % 1 is for the inactive case
        fIndStart = fIndStop + 1;
        fIndStop = fIndStart -1 + numJnodes_i*numStatePJ_i;
        fIndsToLook = fIndStart:numStatePJ_i:fIndStop; % indices of inactive state
        inactiveness_nodes_i = x(fIndsToLook);
        activeStateNodeListInd = find(inactiveness_nodes_i==0);
        if(~isempty(activeStateNodeListInd))
            nodeListInd_i = junctionNodesListInds_i(activeStateNodeListInd);
            nodeIndsActive_i = nodeInds(nodeListInd_i);
            % if any of the active nodes are in the connectionJunction set,
            % make the other nodes in the same set active as well.
            for j=1:numel(nodeIndsActive_i)
                indx = find(connectedJunctionIDs(:,1)==nodeIndsActive_i(j));
                if(~isempty(indx))
                    % this is one of the cluster pixels
                    clusLabel = connectedJunctionIDs(indx,2);
                    clusNodeListInds = find(connectedJunctionIDs(:,2)==clusLabel); 
                    clusNodes = connectedJunctionIDs(clusNodeListInds,1);
                    ilpSegmentation(clusNodes) = 1;
                end
            end
            ilpSegmentation(nodeIndsActive_i) = 1;
        end
    end
end
figure;imagesc(ilpSegmentation);title('ILP contours');
% reconstruct the edges with the values from the orientation filters (HSV)
[output rgbimg] = reconstructHSVgauss_mv(orientedScoreSpace3D,orientations,...
            barLength,barWidth,threshFrac,medianFilterH);
% get the active pixels
output(:,:,3) = ilpSegmentation;
% create HSV image
hsvImage = cat(3,output(:,:,1),output(:,:,2),output(:,:,3));
% convert it to an RGB image
RGBimg = hsv2rgb(hsvImage);
titleStr = sprintf('C = %d : lambda = %d',cNode,decayRate);
figure;imshow(RGBimg);title(titleStr)
