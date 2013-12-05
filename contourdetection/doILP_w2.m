function segmentationOut = doILP_w2(imagePath)

% version 2. 2013.11.29

% using learned param, get segmentation

% ILP script 5
% with the new cost calculation at the junctions, incorporating the
% directionality of the 

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
% imagePath = '/home/thanuja/Dropbox/data/evaldata/input/I11_raw05.tif';

imagePath = '/home/thanuja/Dropbox/data/evaldata/input/I08_raw05.tif';

% labelImagePath = '/home/thanuja/Dropbox/data/evaldata/labels/I05_neuronLabels05.tif';
% imagePath = '/home/thanuja/Dropbox/data/evaldata2/input/I03_raw06.tif';
% labelImagePath = '/home/thanuja/Dropbox/data/evaldata2/labels/I03_neuronLabels06.tif';

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

numTrees = 500;

lenThresh = 25;     % max length of edges to be checked for misorientations
lenThreshBB = 4;    % min length of edges to be considered for being in the backbone (BB)
priorThreshFracBB = 0.55; % threshold of edgePrior for an edge to be considered BB
minNumActEdgesPercentage = 0;  % percentage of the tot num edges to retain (min)

% param
numParam = 0;       % number of params to learn. refer QP section
% edgeOff,edgeOff,nodeoff,nodePos,nodeNeg,regionOff,regionOn

cEdge = 1;        % general scaling factor for edge priors
cCell = 1;        % positive scaling factor for cell priors
cPos = 1;         % scaling factor for positive nodeAngleCosts
cNeg = 1;         % scaling factor for negative nodeAngleCosts
bbEdgeReward = 1;
offEdgeReward = 1;

bbJunctionReward = 1;       % inactivation cost for bbjunction
boundaryEdgeReward = 1;     % prior value for boundary edges so that
                            % they won't have too much weight
                            
% learned parameters
%   1. w_off_e
%   2. w_on_e
%   3. w_off_n
%   4. w_on_n
%   5. w_off_r
%   6. w_on_r

% w_off_e = 1;
% w_on_e = 1;
% w_off_n = 1;
% w_on_n = 1;
% w_off_r = 1;
% w_on_r = 1;

% w_off_e = -11.8404;
% w_on_e = 11.8404;
% w_off_n = 6.33708;
% w_on_n = 10.6658;
% w_off_r = -2.55561;
% w_on_r = 2.55561;

w_off_e = 16.8913;
w_on_e = -16.8913;
w_off_n = -10.6301;
w_on_n = 9.27912;
w_off_r = 1.07527;
w_on_r = -1.07527;



% tot num of int variables = 2*numEdges + 4*numJ3 + 7*numJ4
% coeff (unary prior) for turning off each edge = +edgePriors (col vector)
% coeff (unary prior) for turning on each edge = -edgePriors (col vector)
% coeff for turning off J3s: min(j3NodeAngleCost): max(j3NodeAngleCost)
% coeff for turning on J3-config(1 to 3): j3NodeAngleCost
% coeff for turning off J4s: max(j3NodeAngleCost)
% coeff for turning on J4-config(1 to 6): j4NodeAngleCost

%% read inputimage and get orientedScoreSpace and max_abs value of OFR
disp('using image file:')
disp(imagePath);
imgIn0 = double(imread(imagePath));
% labelImage = imread(labelImagePath);

% add thick border
if(b_imWithBorder)
    imgIn = addThickBorder(imgIn0,marginSize,marginPixVal);
  %   labelImage = addThickBorder(labelImage,marginSize,marginPixVal);
end


%% Oriented Edge Filtering
[output,rgbimg,OFR] = getOFR(imgIn,orientations,...
                        barLength,barWidth,invertImg,threshFrac);
% output is in HSV form
OFR_mag = output(:,:,3);

% generate hsv outputs using the orientation information
% output(:,:,1) contains the hue (orinetation) information
if(~fromInputImage)
    [output rgbimg] = reconstructHSVgauss_mv(OFR,orientations,...
            barLength,barWidth,threshFrac,medianFilterH);
    OFR_mag = imread(imFilePath);
end
if(showIntermediate)
    figure;imshow(rgbimg)
end

%% watershed segmentation
ws = watershed(OFR_mag);
[sizeR,sizeC] = size(ws);
%% generate graph from the watershed edges
disp('creating graph from watershed boundaries...');
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs,selfEdgePixelSet] ...
    = getGraphFromWS(ws,output,showIntermediate);
% nodeEdges - contain edgeIDs for each node

nodeInds = nodeEdges(:,1);                  % indices of the junction nodes
edgeListInds = edges2pixels(:,1);
junctionTypeListInds = getJunctionTypeListInds(nodeEdges);
% col1 has the listInds of type J2, col2: J3 etc. listInds are wrt
% nodeInds list of pixel indices of the detected junctions
clusterNodeIDs = connectedJunctionIDs(:,1); % indices of the clustered junction nodes
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

numEdges = size(edges2nodes,1);

% boundary edges
boundaryEdgeIDs = getBoundaryEdges2(wsRegionBoundariesFromGraph,barLength,edgepixels,...
    nodeEdges,edgeListInds,showIntermediate);
numBoundaryEdges = numel(boundaryEdgeIDs);

[~,boundaryEdgeListInds] = intersect(edgeListInds,boundaryEdgeIDs); 

disp('preparing coefficients for ILP solver...')
%% Edge unary values
% edge priors - from orientation filters
edgePriors = getEdgeUnaryAbs(edgepixels,output(:,:,3));

% get edge activation probabilities from RFC
if ~exist('forestEdgeProb.mat','file')
    disp('RF for edge classification. Training new classifier...')
    forestEdgeProb = trainRF_edgeProb();
else
    load forestEdgeProb.mat
    disp('loaded pre-trained RF for edge activation probability inference.')
end

edgeUnary = getEdgeProbabilitiesFromRFC...
            (forestEdgeProb,imgIn,OFR,edgepixels,edgePriors,...
            boundaryEdgeIDs,edgeListInds,numTrees);


% assigning predetermined edgePriors for boundaryEdges before nodeAngleCost
% calculation
% edgePriors(boundaryEdgeListInds) = boundaryEdgeReward;

% visualize edge unaries
edgeUnaryMat = visualizeEdgeUnaries(edgepixels,edgeUnary,sizeR,sizeC);
if(showIntermediate)
    figure;imagesc(edgeUnaryMat);title('abs-max-OFR')
end
%% Edge pairs - Junction costs
[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);

jEdges = getEdgesForAllNodeTypes(nodeEdges,junctionTypeListInds);
% jEdges{i} - cell array. each cell corresponds to the set of edgeIDs for the
% junction of type i (type1 = J2). A row of a cell corresponds to a node of
% that type of junction.
jAnglesAll = getNodeAnglesForAllJtypes(junctionTypeListInds,...
    nodeInds,jEdges,edges2pixels,OFR,sizeR,sizeC,orientationStepSize);
% jAnglesAll{i} - cell array. each row of a cell corresponds to the set of angles for each
% edge at each junction of type 1 (= J2)

% get the angles for the edges based on its position in the graph
jAnglesAll_alpha = getNodeAngles_fromGraph_allJtypes(junctionTypeListInds,...
    nodeInds,jEdges,edges2pixels,sizeR,sizeC,edges2nodes);

% angle costs
nodeAngleCosts = cell(1,numJtypes);
for i=1:numJtypes
    theta_i = jAnglesAll{i};
    alpha_i = jAnglesAll_alpha{i};
    if(theta_i<0)
        % no such angles for this type of junction
    else
        edgePriors_i = getOrderedEdgePriorsForJ(i,junctionTypeListInds,...
                    nodeEdges,edgeUnary,edgeListInds);
        nodeAngleCosts{i} = getNodeAngleCost_directional(theta_i,alpha_i,...
                                edgePriors_i,w_on_n);
    end
end
%% Faces of wsgraph -> cell types (between pairs of cells)
disp('calculating adjacency graph of regions ...')
[faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] ...
    = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder,boundaryEdgeIDs);

[~,edgeOrientationsInds] = getEdgePriors(OFR,edges2pixels);
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
regionUnary = regionScoreCalculator(forest,normalizedInputImage,setOfRegions,edges2pixels,...
    nodeInds,edges2nodes,cCell,wsIDsForRegions,ws,showIntermediate);
numRegions = numel(regionUnary);
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
%% Get training labels for the image
% 
% [labelImg_indexed,numLabels] = getLabelIndexImg(labelImage);
% [c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,c_extNodeInds]...
%             = getCells2WSregions(labelImg_indexed,ws,numLabels,setOfRegions,...
%             edgeListInds,edges2nodes);
% 
% activeEdgeIDs = getElementsFromCell(c_extEdgeIDs);
% [~,activeEdgeListInds] = intersect(edgeListInds,activeEdgeIDs);
% 
% activeWSregionListInds = getElementsFromCell(c_cells2WSregions);
% 
% activeRegionListInds = activeWSregionListInds - 1;
% 
% activeNodeListInds = getElementsFromCell(c_extNodeInds);
% 
% numEdges = size(edges2nodes,1);

%% visualize training data
% 
% strDataVisualization = visualizeStrData...
%         (c_internalEdgeIDs,c_extEdgeIDs,edgeListInds,edgepixels,...
%         c_internalNodeInds,c_extNodeInds,nodeInds,connectedJunctionIDs,...
%         c_cells2WSregions,ws,numLabels,sizeR,sizeC);
% 
% labelVector = getLabelVector...
%     (activeEdgeListInds,activeNodeListInds,activeRegionListInds,...
%     numEdges,numRegions,jEdges,junctionTypeListInds,edgeListInds);    
% 
% visLV = visualizeXall(labelVector,sizeR,sizeC,numEdges,numRegions,edgepixels,...
%             junctionTypeListInds,nodeInds,connectedJunctionIDs,...
%             nodeEdges,edgeListInds,wsIDsForRegions,ws,twoRegionEdges,edges2regions,...
%             output,showIntermediate);
%         

% TODO:
% labelVectorVisual = visualizeX(labelVector,sizeR,sizeC,numEdges,numRegions,edgepixels,...
%             junctionTypeListInds,nodeInds,connectedJunctionIDs,edges2nodes,...
%             nodeEdges,edgeListInds,faceAdj,setOfRegions,wsIDsForRegions,ws,...
%             marginSize);

%% ILP
% cost function to minimize
% state vector x: {edges*2}{J3*4}{J4*7}

% numJunctions = numel(nodeInds);
% tot num of int variables = 2*numEdges + 4*numJ3 + 7*numJ4
% coeff (unary prior) for turning off each edge = +edgePriors (col vector)
% coeff (unary prior) for turning on each edge = -edgePriors (col vector)
% coeff for turning off J3s: min(j3NodeAngleCost): max(j3NodeAngleCost)
% coeff for turning on J3-config(1 to 3): j3NodeAngleCost
% coeff for turning off J4s: max(j3NodeAngleCost)
% coeff for turning on J4-config(1 to 6): j4NodeAngleCost


% constraints
% equality constraints and closedness constrains in Aeq matrix
% [Aeq,beq,numEq,numLt,numRegionVars] = getConstraints(numEdges,jEdges,edges2pixels,nodeAngleCosts,...
%             offEdgeListIDs,onEdgeListIDs,minNumActEdgesPercentage,...
%             twoRegionEdges,edges2regions,setOfRegions,edgeOrientations,jAnglesAll_alpha,...
%             nodeEdges,junctionTypeListInds,edges2nodes,sizeR,sizeC);
        
[Aeq,beq,senseArray,numEdges,numNodeConf,numRegions,nodeTypeStats]...
    = getILPConstraints(edgeListInds,edges2nodes,nodeEdges,junctionTypeListInds,...
        jEdges,dirEdges2regionsOnOff,setOfRegions);
        

% last 7 variables are continuous RVs corresponding to the parameters to be
% learned.
%   1. w_off_e
%   2. w_on_e
%   3. w_off_n
%   4. w_on_n_neg
%   5. w_on_n_pos
%   6. w_off_r
%   7. w_on_r

% qsparse = getQuadraticObjective_PE(edgePriors,nodeAngleCosts,regionPriors,numParam);
        
% f = getILPcoefficientVector2(scaledEdgePriors,nodeAngleCosts,...
%     bbNodeListInds,junctionTypeListInds,bbJunctionReward,regionPriors);
% numAcols = size(Aeq,2);
% f = zeros(1,numAcols);
% bbJunctionCost = bbJunctionReward;

% f = getILPObjectiveVectorParametric(edgeUnary,nodeAngleCosts,...
%             regionUnary,w_off_e,w_on_e,w_off_n,w_on_n,w_off_r,w_on_r);

f = getILPObjectiveVectorParametric2(edgeUnary_directed,nodeAngleCosts,...
            regionUnary,w_on_e,w_off_n,w_on_n,w_on_r,...
            nodeTypeStats);

% senseArray(1:numEq) = '=';
% if(numLt>0)
%     senseArray((numEq+1):(numEq+numLt)) = '<';
% end
% if(numel(gt_rowID)>0)
%    senseArray(gt_rowID) = '>'; 
% end
% variable types
% vtypeArray(1:numBinaryVar) = 'B'; % binary
% vtypeArray((numBinaryVar+1):(numBinaryVar+numParam)) = 'C'; % continuous
% lower bounds
% lbArray(1:(numBinaryVar+numParam)) = 0;
% upper bounds
% ubArray(1:(numBinaryVar+numParam)) = 1;

%% Write files for structured learninig bmrm
% featureMat = writeFeaturesFile(f,jEdges,numEdges,numRegions);
% 
% constraints = writeConstraintsFile(Aeq,beq,senseArray);
% 
% features = writeLabelsFile(labelVector);

%% solver
if(useGurobi)
    disp('using Gurobi ILP solver...');
    model.A = sparse(Aeq);
    model.rhs = beq;
    model.obj = f';
    model.sense = senseArray;
    % model.vtype = vtypeArray;
    model.vtype = 'B';
    % model.lb = lbArray;
    % model.ub = ubArray;
    model.modelname = 'contourDetectionILP1';
    % initial guess
    % model.start = labelVector;
    
    
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