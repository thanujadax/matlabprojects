% function segmentationOut = doILPseg(imagePath)

% ILP script 5
% with the new cost calculation at the junctions, incorporating the
% directionality of the 

gridResolution = 6;

showIntermediate = 0;
useGurobi = 1;
fromInputImage = 1;
% imagePath = '/home/thanuja/Dropbox/data/mitoData/emJ_00_170x.png';
% imagePath = '/home/thanuja/Dropbox/data/testImg/testCurves1.png';
% imagePath = '/home/thanuja/Dropbox/data/mitoData/stem1_256by256.png';
% imagePath = '/home/thanuja/Dropbox/data/thanuja/emchallenge-class/competition-final0000.tif';
% hard coded back bone edge 1962
% imagePath = '/home/thanuja/Dropbox/data/RF_training_edge/I15_testingImage.tif';
% imagePath = '/home/thanuja/Dropbox/data/RF_training_edge/I05_trainingImage.tif';
imagePath = '/home/thanuja/Dropbox/data/evaldata/input/I10_raw05.tif';
% imagePath = '/home/thanuja/Dropbox/data/edgeTraining2/trainingRaw/I00_raw05.tif';

orientationStepSize = 10;
orientations = 0:orientationStepSize:350;

barLength = 13; % should be odd
barWidth = 4; %
marginSize = ceil(barLength/2);
marginPixVal = 0;
threshFrac = 0.1;
medianFilterH = 0;
invertImg = 1;      % 1 for EM images when input image is taken from imagePath
b_imWithBorder = 1; % add thick dark border around the image
eps_orientation = 25; % for edge response calculation from OFR.
                    % range of orientations = orientation +- eps_orientation
barHalfWidth = 2;  % for edge response calculation from OFR. Additional thickness
                % on each side of the edge to consider OFR around it.
% param
cEdge = 10;        % general scaling factor for edge priors
% cNode = 100;        % scaling factor for the node cost coming from gaussian normal distr.
cCell = 100;        % positive scaling factor for cell priors
% sig = 50;         % standard deviation(degrees) for the node cost function's gaussian distr.
% midPoint = 180;   % angle difference of an edge pair (in degrees) for maximum cost 
lenThresh = 25;     % max length of edges to be checked for misorientations
lenThreshBB = 4;    % min length of edges to be considered for being in the backbone (BB)
priorThreshFracBB = 0.55; % threshold of edgePrior for an edge to be considered BB

% cPos = 1000;        % scaling factor for positive nodeAngleCosts
% cNeg = 10;          % scaling factor for negative nodeAngleCosts
w_on_n = 100;       % scaling factor for nodeAngleCosts

minNumActEdgesPercentage = 0;  % percentage of the tot num edges to retain (min)
bbEdgeReward = 1500;
offEdgeReward = -500;
bbJunctionReward = 1000;        % inactivation cost for bbjunction
boundaryEdgeReward = -35;   % prior value for boundary edges so that they won't have too much weight


%% read inputimage and get orientedScoreSpace and max_abs value of OFR
disp('using image file:')
disp(imagePath);
imgIn0 = double(imread(imagePath));

% add border
if(b_imWithBorder)
    imgIn = addThickBorder(imgIn0,marginSize,marginPixVal);
end


%% Oriented Edge Filtering
[output,rgbimg,orientedScoreSpace3D] = getOFR(imgIn,orientations,...
                        barLength,barWidth,invertImg,threshFrac);
% output is in HSV form
OFR_mag = output(:,:,3);

% generate hsv outputs using the orientation information
% output(:,:,1) contains the hue (orinetation) information
if(~fromInputImage)
    [output rgbimg] = reconstructHSVgauss_mv(orientedScoreSpace3D,orientations,...
            barLength,barWidth,threshFrac,medianFilterH);
    OFR_mag = imread(imFilePath);
end
if(showIntermediate)
    figure;imshow(rgbimg)
end

%% watershed segmentation
%ws = watershed(OFR_mag);
[ws,setOfRegions,edges2pixels,edges2nodes,nodeEdges,adjacencyMat,nodeInds,...
    edges2regions,boundaryEdgeListInds,twoRegionEdges,faceAdj,wsIDsForRegions]...
                            = getImageGrid(imgIn,gridResolution,showIntermediate);
                                    
[sizeR,sizeC] = size(ws);
connectedJunctionIDs = [];
selfEdgePixelSet = [];

numRegions = size(setOfRegions,1);

%% generate graph from the watershed edges
disp('creating graph from watershed boundaries...');
% [adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs,selfEdgePixelSet] ...
%     = getGraphFromWS(ws,output,showIntermediate);
% nodeInds = nodeEdges(:,1);                  % indices of the junction nodes

edgeListInds = edges2pixels(:,1);
junctionTypeListInds = getJunctionTypeListInds(nodeEdges);

% col1 has the listInds of type J2, col2: J3 etc. listInds are wrt
% nodeInds list of pixel indices of the detected junctions
% clusterNodeIDs = connectedJunctionIDs(:,1); % indices of the clustered junction nodes
clusterNodeIDs = [];

disp('graph created!')
wsRegionBoundariesFromGraph = zeros(sizeR,sizeC);
wsRegionBoundariesFromGraph(nodeInds) = 0.7;          % junction nodes
wsRegionBoundariesFromGraph(clusterNodeIDs) = 0.5;    % cluster nodes
[nre,nce] = size(edges2pixels);  % first column is the edgeID
edgepixels = edges2pixels(:,2:nce);
wsRegionBoundariesFromGraph(edgepixels(edgepixels>0)) = 1; % edge pixels
if(showIntermediate)
    figure;imagesc(wsRegionBoundariesFromGraph);title('boundaries from graph') 
end
% boundary edges
% boundaryEdgeIDs = getBoundaryEdges2(wsRegionBoundariesFromGraph,barLength,edgepixels,...
%     nodeEdges,edgeListInds,showIntermediate);
numBoundaryEdges = numel(boundaryEdgeListInds);

% [~,boundaryEdgeListInds] = intersect(edgeListInds,boundaryEdgeListInds); 


disp('preparing coefficients for ILP solver...')
%% Edge unary values
jEdges = getEdgesForAllNodeTypes(nodeEdges,junctionTypeListInds);
% jEdges{i} - cell array. each cell corresponds to the set of edges for the
% junction of type i (type1 = J2). A row of a cell corresponds to a node of
% that type of junction.

% jAnglesAll = getNodeAnglesForAllJtypes(junctionTypeListInds,...
%     nodeInds,jEdges,edges2pixels,orientedScoreSpace3D,sizeR,sizeC,orientationStepSize);

% jAnglesAll{i} - cell array. each row of a cell corresponds to the set of angles for each
% edge at each junction of type 1 (= J2)

[edgeOrientationsList, edgeResponses] = getEdgeOrientations_grid...
            (edges2nodes,sizeR,sizeC,eps_orientation,orientedScoreSpace3D,...
            edges2pixels,gridResolution,barHalfWidth,orientationStepSize,...
            nodeInds);
% edgeOrientationsList: col1- default orientation, col2-complementary
% orientation
% edgeResponses: edge response (signed) to the default orientation

jAnglesAll = getNodeAnglesFromEdgeOrientations...
            (jEdges,edgeOrientationsList);

% get the angles for the edges based on its position in the graph
jAnglesAll_alpha = getNodeAngles_fromGraph_allJtypes(junctionTypeListInds,...
    nodeInds,jEdges,edges2pixels,sizeR,sizeC,edges2nodes);

% edge priors - from orientation filters
% edgePriors = getEdgeUnaryAbs(edgepixels,output(:,:,3));
edgePriors = edgeResponses .* 10;

% assigning predetermined edgePriors for boundaryEdges before nodeAngleCost
% calculation
% edgePriors(boundaryEdgeListInds) = boundaryEdgeReward;

% visualize edge unaries
edgeUnaryMat = visualizeEdgeUnaries(edgepixels,edgePriors,sizeR,sizeC);
if(showIntermediate)
    figure;imagesc(edgeUnaryMat);title('abs-max-OFR')
end
%% Edge pairs - Junction costs
[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);

% angle costs
nodeAngleCosts = cell(1,numJtypes);
for i=1:numJtypes
    theta_i = jAnglesAll{i};
    alpha_i = jAnglesAll_alpha{i};
    if(theta_i<0)
        % no such angles for this type of junction
    else
        edgePriors_i = getOrderedEdgePriorsForJ(i,junctionTypeListInds,...
                    nodeEdges,edgePriors,edgeListInds);
        nodeAngleCosts{i} = getNodeAngleCost_directional(theta_i,alpha_i,...
                                edgePriors_i,w_on_n);
    end
end
%% Faces of wsgraph -> cell types (between pairs of cells)
disp('calculating adjacency graph of regions ...')
% [faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] ...
%     = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder,boundaryEdgeListInds);

[~,edgeOrientationsInds] = getEdgePriors(orientedScoreSpace3D,edges2pixels);
edgeOrientations = (edgeOrientationsInds-1).*orientationStepSize;

% normalize input image
normalizedInputImage = imgIn./(max(max(imgIn)));

% get cell priors from RFC probability map
if ~exist('forest.mat','file')
    disp('RF for membrane classification not found. Training new classifier...')
    forest = trainRandomForest_pixelProb();
else
    load forest.mat
    disp('loaded pre-trained RF for membrane vs cell-interior classification')
end
regionPriors = regionScoreCalculator(forest,normalizedInputImage,setOfRegions,edges2pixels,...
    nodeInds,edges2nodes,cCell,wsIDsForRegions,ws,showIntermediate);
%% Boundary edges
% assigning predetermined edge priors for boundary edges after
% nodeAngleCost calculation

edgePriors(boundaryEdgeListInds) = boundaryEdgeReward;

% boundaryNodeEdges = all edges that has at least one end connected to the boundary
boundaryNodeListInds = edges2nodes(boundaryEdgeListInds,:);
boundaryNodeListInds = unique(boundaryNodeListInds);
boundaryNodeEdges = nodeEdges(boundaryNodeListInds,:);
boundaryNodeEdges(:,1)=[];
boundaryNodeEdges = unique(boundaryNodeEdges);
boundaryNodeEdges = boundaryNodeEdges(boundaryNodeEdges>0);
boundaryNodeEdgeListIDs = numel(boundaryNodeEdges);
for i=1:numel(boundaryNodeEdges)
    boundaryNodeEdgeListIDs(i) = find(edgeListInds==boundaryNodeEdges(i));
end

%% Removing misoriented edges
% uses the compatibility of the orientation of the adjoining pixels of each
% edge

offEdgeListIDs = getUnOrientedEdgeIDs(edgepixels,...
                lenThresh,output(:,:,1),sizeR,sizeC);
            
% remove boundaryNodeEdgeListIDs from the offEdgeListIDs
offEdgeListIDs = setdiff(offEdgeListIDs,boundaryNodeEdgeListIDs);
          
% setting edgePriors
edgePriors(offEdgeListIDs) = offEdgeReward;
% visualize off edges
imgOffEdges = visualizeOffEdges(offEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC);
if(showIntermediate)
    figure;imshow(imgOffEdges); title('visualization of edges turned off')
end
%% Backbone

onEdgeListIDs = getBackboneEdgeIDs(edgepixels,edgePriors,...
                lenThreshBB,priorThreshFracBB);
% removing boundary edges from the onEdgeListID list
onEdgeListIDs = setdiff(onEdgeListIDs,boundaryEdgeListInds);

bbNodeListInds = getJunctionsForEdges(edges2nodes,onEdgeListIDs);
bbNodePixInds = nodeInds(bbNodeListInds);
bbnodesVis = wsRegionBoundariesFromGraph;
bbnodesVis(bbNodePixInds) = 0.2;
if(showIntermediate)
    figure;imagesc(bbnodesVis);
end
% visualize BB edges
imgBBEdges = visualizeOffEdges(onEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC);
if(showIntermediate)
    figure;imshow(imgBBEdges); title('visualization of backbone 1')
end
edgePriors(onEdgeListIDs) = bbEdgeReward;
% test - enforcing some edges to be picked
clear onEdgeListIDs
onEdgeListIDs = [];

% visualize BB edges
imgBBEdges = visualizeOffEdges(onEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC);
if(showIntermediate)
    figure;imshow(imgBBEdges); title('visualization of backbone 2')
end
    % onEdgeListIDs = 2024;
% % visualize BB edges
% imgBBEdges = visualizeOffEdges(onEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC);
% figure;imshow(imgBBEdges); title('visualization of backbone - constr - 2 ')
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

% constraints
% equality constraints and closedness constrains in Aeq matrix
[Aeq,beq,numEq,numLt,numRegionVars] = getConstraints(numEdges,jEdges,edges2pixels,nodeAngleCosts,...
            offEdgeListIDs,onEdgeListIDs,minNumActEdgesPercentage,...
            twoRegionEdges,edges2regions,setOfRegions,edgeOrientations,jAnglesAll_alpha,...
            nodeEdges,junctionTypeListInds,edges2nodes,sizeR,sizeC);
        
f = getILPcoefficientVector2(scaledEdgePriors,nodeAngleCosts,...
    bbNodeListInds,junctionTypeListInds,bbJunctionReward,regionPriors);
                
senseArray(1:numEq) = '=';
if(numLt>0)
    senseArray((numEq+1):(numEq+numLt)) = '<';
end
%% solver
if(useGurobi)
    disp('using Gurobi ILP solver...');
    model.A = sparse(Aeq);
    model.rhs = beq;
    model.obj = f';
%     model.sense = '=';  % for the constraints given in A
    model.sense = senseArray;
    model.vtype = 'B';  % binary variables. otherwise we can specify variable 
    % types for each column of A. e.g. continuous, integer etc
    model.modelname = 'contourDetectionILP1';
    
    params.LogFile = 'gurobi.log';
    params.Presolve = 0;
    params.ResultFile = 'modelfile.mps';
    params.InfUnbdInfo = 1;

    resultGurobi = gurobi(model,params);
    x = resultGurobi.x;
    
    
else
    % Matlab ILP solver
    disp('using MATLAB ILP solver...');
    Initial values for the state variables
    x0 = getInitValues(numEdges,numJ3,numJ4);  % TODO: infeasible. fix it!!
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
segmentationOut = visualizeX(x,sizeR,sizeC,numEdges,numRegions,edgepixels,...
            junctionTypeListInds,nodeInds,connectedJunctionIDs,edges2nodes,...
            nodeEdges,edgeListInds,faceAdj,setOfRegions,wsIDsForRegions,ws,...
            marginSize,showIntermediate);