function segmentationOut = doILP_w_dir(inputPath,imageFileName,imageID,...
    rawType,neuronProbabilityType,membraneProbabilityType,mitoProbabilityType)

% version 3. 2014.01.06

% each edge in the ws graph is represented by 2 (oppositely) directed edges 

%% Settings

produceBMRMfiles = 0;  % set label file below if 1
showIntermediate = 1;
fillSpaces = 1;          % fills holes in segmentationOut
useGurobi = 1;
fromInputImage = 1;
usePrecomputedProbabilityMaps = 1;
useMitochondriaDetection = 0;

%% File names and paths

% inputPath = '/home/thanuja/Dropbox/data2';

% probability map image file should have the same name as the raw image file
rawImagePath = fullfile(inputPath,'raw');
% rawImagePath = '/home/thanuja/Dropbox/data2/raw';
% rawImagePath = '/home/thanuja/Dropbox/data2/probabilities/neuron';
% rawImageFileName = '00.png';
if(~isempty(imageID))
    imgFileString = strcat('*.',rawType);
    rawImageFilesAll = dir(fullfile(rawImagePath,imgFileString));
    rawImageFileName = rawImageFilesAll(imageID).name;

else
    rawImageFileName = imageFileName;
end

rawImageFullFile = fullfile(rawImagePath,rawImageFileName);

% probabilityMapPath = '/home/thanuja/Dropbox/data2/probabilities';
% probabilityMapPath = fullfile(inputPath,'probabilities');
probabilityMapPath = inputPath;

dir_membraneProb = 'membranes';
dir_mitochondriaProb = 'mitochondria';
dir_neuronProb = 'neurons';


if(~isempty(imageID))
    imgFileString = strcat('*.',membraneProbabilityType);
    
    membraneProbabilityImageFilesAll = dir(fullfile(...
        probabilityMapPath,dir_membraneProb,imgFileString));
    membraneProbabilityImage = membraneProbabilityImageFilesAll(imageID).name;
  
    imgFileString = strcat('*.',neuronProbabilityType);
    neuronProbabilityImageFilesAll = dir(fullfile(...
        probabilityMapPath,dir_neuronProb,imgFileString));
    neuronProbabilityImage = neuronProbabilityImageFilesAll(imageID).name;
    neuronProbabilityImage = fullfile(probabilityMapPath,dir_neuronProb,neuronProbabilityImage);
    disp(neuronProbabilityImage);
    
    if(useMitochondriaDetection)
    imgFileString = strcat('*.',mitoProbabilityType);
    mitoProbabilityImageFilesAll = dir(fullfile(...
        probabilityMapPath,dir_mitochondriaProb,imgFileString));
    mitochondriaProbabilityImage = mitoProbabilityImageFilesAll(imageID).name;
    mitochondriaProbabilityImage = fullfile(...
        probabilityMapPath,dir_mitochondriaProb,mitochondriaProbabilityImage);
    disp(mitochondriaProbabilityImage)
    else
        mitochondriaProbabilityImage = [];
    end
    
else
    membraneProbabilityImage = fullfile(probabilityMapPath,dir_membraneProb,rawImageFileName);
    neuronProbabilityImage = fullfile(probabilityMapPath,dir_neuronProb,rawImageFileName);
    if(useMitochondriaDetection)
        mitochondriaProbabilityImage = fullfile(probabilityMapPath,dir_mitochondriaProb,rawImageFileName);
    end
end
% for sbmrm
if(produceBMRMfiles)
    labelImagePath = '/home/thanuja/Dropbox/data2/results';
    labelImageFileName = '00.png';
    labelImageFullFile = fullfile(labelImagePath,labelImageFileName);
end

%% Parameters
orientationStepSize = 10;
orientations = 0:orientationStepSize:350;

barLength = 13; % should be odd
barWidth = 4; %
marginSize = ceil(barLength/2);
marginPixVal = 0.3;
threshFrac = 0.1;   % 0.1 for raw images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
medianFilterH = 0;
invertImg = 1;      % 1 for EM images when input image is taken from imagePath
b_imWithBorder = 1; % add thick dark border around the image

numTrees = 500;

lenThresh = 25;     % max length of edges to be checked for misorientations
lenThreshBB = 4;    % min length of edges to be considered for being in the backbone (BB)
priorThreshFracBB = 0.55; % threshold of edgePrior for an edge to be considered BB
minNumActEdgesPercentage = 0;  % percentage of the tot num edges to retain (min)

% param

cEdge = 1;        % general scaling factor for edge priors
cCell = 1;        % positive scaling factor for cell priors
cPos = 1;         % scaling factor for positive nodeAngleCosts
cNeg = 1;         % scaling factor for negative nodeAngleCosts
bbEdgeReward = 1;
offEdgeReward = 1;

bbJunctionReward = 1;       % inactivation cost for bbjunction
boundaryEdgeReward = 1;     % prior value for boundary edges so that
                            % they won't have too much weight
%% Set parameters
regionOffThreshold = 0.21;  % threshold to pick likely off regions to set off score 
% to +1*w_off_r in the objective
                            
%%  learned parameters
%   1. w_off_e
%   2. w_on_e
%   3. w_off_n
%   4. w_on_n
%   5. w_off_r
%   6. w_on_r
if(produceBMRMfiles)
    % set all parameters to  be learned to 1
    w_on_e = 1;     % edge weight
    w_off_e = 1;
    w_off_n = 1;    % node off weight
    w_on_n = 1;     % node on weight
    w_on_r = 1;     % region weight
    w_off_r = 1;
else
    % use pre-learned parameters
    % optimial w is [-7.52064, -7.38296, 0.468054, 0.403942, -7.79221, -5.75401]
%    w_on_e = -7.52064;     % edge weight
    w_on_e = -10;
    w_off_e = -7.38296;

    w_off_n = 0.468054;    % node off weight
    w_on_n = 0.403942;     % node on weight
    w_on_r = -7.79221;     % region weight
    % w_off_r = -5.75401;
    w_off_r = -18;
end


% tot num of int variables = 2*numEdges + 4*numJ3 + 7*numJ4
% coeff (unary prior) for turning off each edge = +edgePriors (col vector)
% coeff (unary prior) for turning on each edge = -edgePriors (col vector)
% coeff for turning off J3s: min(j3NodeAngleCost): max(j3NodeAngleCost)
% coeff for turning on J3-config(1 to 3): j3NodeAngleCost
% coeff for turning off J4s: max(j3NodeAngleCost)
% coeff for turning on J4-config(1 to 6): j4NodeAngleCost

%% read inputimage and get orientedScoreSpace and max_abs value of OFR
disp('using image file:')
disp(rawImageFullFile);
imgIn0 = double(imread(rawImageFullFile));
[a,b,c] = size(imgIn0);
if(c==3)
    imgIn0 = rgb2gray(imgIn0);
end

% imgIn0 = imgIn0(1:128,:);

if(produceBMRMfiles)
    labelImage = imread(labelImageFullFile);
    % labelImage = labelImage(1:128,:,:);
end
% add thick border
if(b_imWithBorder)
    imgIn = addThickBorder(imgIn0,marginSize,marginPixVal);
    if(produceBMRMfiles)
        labelImage = addThickBorder(labelImage,marginSize,marginPixVal);
    end
end


%% Oriented Edge Filtering
[output,rgbimg,OFR] = getOFR(imgIn,orientations,...
                        barLength,barWidth,invertImg,threshFrac);
% output is in HSV form
OFR_mag = output(:,:,3);
OFR_hue = output(:,:,1);
% generate hsv outputs using the orientation information
% output(:,:,1) contains the hue (orinetation) information

if(showIntermediate)
    figure;imshow(rgbimg)
end

%% watershed segmentation
ws = watershed(OFR_mag);
[sizeR,sizeC] = size(ws);
% randomize WS region IDs
ws = assignRandomIndicesToWatershedTransform(ws);
%% generate graph from the watershed edges
disp('creating graph from watershed boundaries...');
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs,selfEdgePixelSet,...
    ws,ws_original,removedWsIDs,newRemovedEdgeLIDs] ...
    = getGraphFromWS(ws,output,showIntermediate);
clear adjacencyMat
clear output
% nodeEdges - contain edgeIDs for each node

nodeInds = nodeEdges(:,1);                  % indices of the junction nodes
edgeListInds = edges2pixels(:,1);
junctionTypeListInds = getJunctionTypeListInds(nodeEdges);
% col1 has the listInds of type J2, col2: J3 etc. listInds are wrt
% nodeInds list of pixel indices of the detected junctions
if(size(connectedJunctionIDs,2)==2)
    clusterNodeIDs = connectedJunctionIDs(:,1); % indices of the clustered junction nodes
else
    clusterNodeIDs = 0;
end
disp('graph created!')
wsRegionBoundariesFromGraph = zeros(sizeR,sizeC);
wsRegionBoundariesFromGraph(nodeInds) = 0.7;          % junction nodes
if(size(connectedJunctionIDs,2)==2)
    wsRegionBoundariesFromGraph(clusterNodeIDs) = 0.5;    % cluster nodes
end
[nre,nce] = size(edges2pixels);  % first column is the edgeID
edgepixels = edges2pixels(:,2:nce);
wsRegionBoundariesFromGraph(edgepixels(edgepixels>0)) = 1; % edge pixels
if(showIntermediate)
    figure;imagesc(wsRegionBoundariesFromGraph);title('boundaries from graph') 
end

numEdges = size(edges2nodes,1);

% boundary edges
% boundaryEdgeIDs = getBoundaryEdges2(wsRegionBoundariesFromGraph,barLength,edgepixels,...
%     nodeEdges,edgeListInds,showIntermediate);

boundaryEdgeIDs = getBoundaryEdgeIDs(ws,edges2pixels);
numBoundaryEdges = numel(boundaryEdgeIDs);

[~,boundaryEdgeListInds] = intersect(edgeListInds,boundaryEdgeIDs); 

disp('preparing coefficients for ILP solver...')
%% Edge unary values
% edge priors - from orientation filters
edgePriors = getEdgeUnaryAbs(edgepixels,OFR_mag);

% get edge activation probabilities from RFC

if(0) % not using precomputed probability maps for graph edges - doesn't make sense!
    % calculate edgeUnary from probability map image
    edgeUnary = getEdgeProbabilityFromMap(...
        membraneProbabilityImage,edgepixels,marginSize,(1-marginPixVal));
else
    
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
end

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
    nodeInds,jEdges,edges2pixels,sizeR,sizeC,edges2nodes,connectedJunctionIDs);

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
%% Faces of wsgraph -> region types (between pairs of regions)
disp('calculating adjacency graph of regions ...')
[faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] ...
    = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder,boundaryEdgeIDs);

[~,edgeOrientationsInds] = getEdgePriors(OFR,edges2pixels);
clear OFR
edgeOrientations = (edgeOrientationsInds-1).*orientationStepSize;

% normalize input image
normalizedInputImage = imgIn./(max(max(imgIn)));

%% get region unaries
if(usePrecomputedProbabilityMaps)
    
    regionUnary = getRegionScoreFromProbImage(...
    neuronProbabilityImage,mitochondriaProbabilityImage,...
    useMitochondriaDetection,marginSize,marginPixVal,...
    setOfRegions,sizeR,sizeC,wsIDsForRegions,ws,showIntermediate);
    
else
    
    if ~exist('forest.mat','file')
        disp('RF for membrane classification not found. Training new classifier...')
        forest = trainRandomForest_pixelProb();
    else
        load forest.mat
        disp('loaded pre-trained RF for membrane vs cell-interior classification')
    end
    regionUnary = regionScoreCalculator(forest,normalizedInputImage,setOfRegions,edges2pixels,...
        nodeInds,edges2nodes,cCell,wsIDsForRegions,ws,showIntermediate);   
    
end


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
                lenThresh,OFR_hue,sizeR,sizeC);
            
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

% commenting out the following two methods : 20150518
% [c_edgeLIDsForRegions_dir_cw,setOfRegions_edgeLIDs,edgeLIDs2nodes_directional] ...
%         = getOrderedRegionEdgeListIndsDir...
%         (setOfRegions,edges2nodes,jAnglesAll_alpha,...
%         junctionTypeListInds,nodeEdges,edgeListInds,edges2pixels,sizeR,sizeC);
% 
% dirEdges2regionsOnOff = getRegionsForDirectedEdges...
%             (c_edgeLIDsForRegions_dir_cw,edgeLIDs2nodes_directional,...
%             setOfRegions_edgeLIDs,numEdges);
% dEdges2regionsOnOff = edgeListInd_dir (=rowID) | onRegion | offRegion  : dir N1->N2
%   regionID = 0 is for the image border.

% TODO : unfinished implementation 20150518
dirEdges2regionsOnOff ...
        = getDirectedEdgeOnOffRegions...
        (setOfRegions,edges2nodes,jAnglesAll_alpha,...
        junctionTypeListInds,nodeEdgeIDs,edgeListIndsAll,...
        edges2pixels,sizeR,sizeC);

if(produceBMRMfiles)
    [labelImg_indexed,numLabels] = getLabelIndexImg(labelImage);
    [c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,c_extNodeInds]...
                = getCells2WSregions(labelImg_indexed,ws,numLabels,setOfRegions,...
                edgeListInds,edges2nodes);
    activeWSregionListInds_tr = getElementsFromCell(c_cells2WSregions);  
else
    activeWSregionListInds_tr = [];
end

[model.A,b,senseArray,numEdges,numNodeConf,numRegions,nodeTypeStats]...
    = getILPConstraints(edgeListInds,edges2nodes,nodeEdges,junctionTypeListInds,...
        jEdges,dirEdges2regionsOnOff,setOfRegions,activeWSregionListInds_tr);
        

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

% f = getILPObjectiveVectorParametric2(edgeUnary_directed,nodeAngleCosts,...
%             regionUnary,w_on_e,w_off_n,w_on_n,w_on_r,...
%             nodeTypeStats);
        
if(produceBMRMfiles)
    f = getILPObjVect_Tr(labelImage,ws,edgeListInds,...
                setOfRegions,edges2nodes,numEdges,numNodeConf,numRegions,...
                edgeUnary);
else            
    f = getILPObjectiveVectorParametric2(edgeUnary,nodeAngleCosts,...
            regionUnary,w_on_e,w_off_e,w_off_n,w_on_n,w_on_r,w_off_r,...
            nodeTypeStats,offEdgeListIDs,regionOffThreshold,numNodeConf);
end

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


%% solver
if(useGurobi)
    disp('using Gurobi ILP solver...');
    % model.A = sparse(double(A));
    model.rhs = b;
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
%% write BMRM files
if (produceBMRMfiles)
    f = getILPObjectiveVectorParametric2(edgeUnary,nodeAngleCosts,...
            regionUnary,w_on_e,w_off_e,w_off_n,w_on_n,w_on_r,w_off_r,...
            nodeTypeStats,offEdgeListIDs,regionOffThreshold,numNodeConf); % w's are set to 1.
    featureMat = writeFeaturesFile2(f,jEdges,numEdges,numRegions);
    constraints = writeConstraintsFile(model.A,b,senseArray);
    labels = writeLabelsFile(x);
end

%% visualize
segmentationOut = visualizeX2(x,sizeR,sizeC,numEdges,numRegions,edgepixels,...
            junctionTypeListInds,nodeInds,connectedJunctionIDs,edges2nodes,...
            nodeEdges,edgeListInds,faceAdj,setOfRegions,wsIDsForRegions,ws,...
            marginSize,showIntermediate,fillSpaces);