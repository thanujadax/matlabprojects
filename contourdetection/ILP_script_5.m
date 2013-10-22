% ILP script 5
% with the new cost calculation at the junctions, incorporating the
% directionality of the 

isToyProb = 0;
useGurobi = 1;
fromInputImage = 1;
% imagePath = '/home/thanuja/Dropbox/data/mitoData/emJ_00_170x.png';
% imagePath = '/home/thanuja/Dropbox/data/testImg/testCurves1.png';
% imagePath = '/home/thanuja/Dropbox/data/mitoData/stem1_256by256.png';
% imagePath = '/home/thanuja/Dropbox/data/thanuja/emchallenge-class/competition-final0000.tif';
imagePath = '/home/thanuja/Dropbox/data/RF_training_edge/I15_testingImage.tif';
% imagePath = '/home/thanuja/Dropbox/data/RF_training_edge/I05_trainingImage.tif';
orientationsStepSize = 10;
orientations = 0:orientationsStepSize:350;

barLength = 13; % should be odd
barWidth = 4; %
marginSize = ceil(barLength/2);
% marginPixVal = 0.1;
addBorder = ceil(barLength/2);
threshFrac = 0.1;
medianFilterH = 0;
invertImg = 1;      % 1 for EM images when input image is taken from imagePath
% max vote response image of the orientation filters
b_imWithBorder = 1; % add thick dark border around the image

if(isToyProb)
%     imFilePath = 'testMem4_V.png';
    imFilePath = 'circle1_V.png';
    % votes for each orientation for each edge
%     load('orientedScoreSpace3D.mat') % loads the orientation filter scores
    load('orientedScoreSpace3D_circle1.mat') % loads the orientation filter scores
elseif(fromInputImage)
    % read inputimage and get orientedScoreSpace and max_abs value of OFR
    disp('using image file:')
    disp(imagePath);
    imgIn0 = double(imread(imagePath));
    % add border
    marginPixVal = min(min(imgIn0));
    if(b_imWithBorder)
        imgIn = addThickBorder(imgIn0,marginSize,marginPixVal);
    end
    [output,rgbimg,orientedScoreSpace3D] = getOFR(imgIn,orientations,...
                            barLength,barWidth,invertImg,threshFrac);
    % output is in HSV form
    OFR_mag = output(:,:,3);
else
%     imFilePath = 'stem_256x_t02_V.png';
    imFilePath = '/home/thanuja/Dropbox/data/mitoData/emJ_00_350x_V.png';
    % votes for each orientation for each edge
%     load('orientedScoreSpace3D_stem256x.mat') % loads the orientation filter scores
    load('orientedScoreSpace3D_emJ350x.mat')
end

angleStep = 10; % 10 degrees discretization step of orientations

% param
cEdge = 1;        % general scaling factor for edge priors
% cNode = 100;        % scaling factor for the node cost coming from gaussian normal distr.
cCell = 1000;        % positive scaling factor for cell priors
% sig = 50;         % standard deviation(degrees) for the node cost function's gaussian distr.
% midPoint = 180;   % angle difference of an edge pair (in degrees) for maximum cost 
lenThresh = 25;     % max length of edges to be checked for misorientations
lenThreshBB = 4;    % min length of edges to be considered for being in the backbone (BB)
priorThreshFracBB = 0.55; % threshold of edgePrior for an edge to be considered BB

maxCost_direction = 1000;  % C for the directional cost function
cPos = 1000;        % scaling factor for positive nodeAngleCosts
cNeg = 10;          % scaling factor for negative nodeAngleCosts
minNumActEdgesPercentage = 0;  % percentage of the tot num edges to retain (min)
bbEdgeReward = 1500;
offEdgeReward = -500;
bbJunctionReward = 1000;        % inactivation cost for bbjunction
boundaryEdgeReward = -35;   % prior value for boundary edges so that they won't have too much weight


% generate hsv outputs using the orientation information
% output(:,:,1) contains the hue (orinetation) information
if(~fromInputImage)
    [output rgbimg] = reconstructHSVgauss_mv(orientedScoreSpace3D,orientations,...
            barLength,barWidth,threshFrac,medianFilterH);
    OFR_mag = imread(imFilePath);
end
figure;imshow(rgbimg)
OFR_abs = output(:,:,3);
% watershed segmentation
ws = watershed(OFR_mag);
% ws = watershed(output(:,:,3));
[sizeR,sizeC] = size(ws);
%% generate graph from the watershed edges
disp('creating graph from watershed boundaries...');
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs,selfEdgePixelSet] ...
    = getGraphFromWS(ws,output);
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
figure;imagesc(wsRegionBoundariesFromGraph);title('boundaries from graph') 

% boundary edges
boundaryEdges = getBoundaryEdges2(wsRegionBoundariesFromGraph,barLength,edgepixels,...
    nodeEdges,edgeListInds);
numBoundaryEdges = numel(boundaryEdges);
boundaryEdgeListInds = zeros(numBoundaryEdges,1);
for i=1:numBoundaryEdges
    boundaryEdgeListInds(i) = find(edgeListInds==boundaryEdges(i));
end

disp('preparing coefficients for ILP solver...')
%% Edge unary values
% edge priors - from orientation filters
% edgePriors = getEdgePriors(orientedScoreSpace3D,edges2pixels);
edgePriors = getEdgeUnaryAbs(edgepixels,output(:,:,3));

% assigning predetermined edgePriors for boundaryEdges before nodeAngleCost
% calculation
% edgePriors(boundaryEdgeListInds) = boundaryEdgeReward;

% visualize edge unaries
edgeUnaryMat = visualizeEdgeUnaries(edgepixels,edgePriors,sizeR,sizeC);
figure;imagesc(edgeUnaryMat);title('abs-max-OFR')
%% Edge pairs - Junction costs
[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);

jEdges = getEdgesForAllNodeTypes(nodeEdges,junctionTypeListInds);
% jEdges{i} - cell array. each cell corresponds to the set of edges for the
% junction of type i (type1 = J2). A row of a cell corresponds to a node of
% that type of junction.
jAnglesAll = getNodeAnglesForAllJtypes(junctionTypeListInds,...
    nodeInds,jEdges,edges2pixels,orientedScoreSpace3D,sizeR,sizeC,angleStep);
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
        %nodeAngleCosts{i} = getNodeAngleCost(dTheta_i,midPoint,sig,cNode);
        edgePriors_i = getOrderedEdgePriorsForJ(i,junctionTypeListInds,...
                    nodeEdges,edgePriors,edgeListInds);
        nodeAngleCosts{i} = getNodeAngleCost_directional(theta_i,alpha_i,...
                                edgePriors_i,cPos,cNeg);
    end
end
%% Faces of wsgraph -> cell types (between pairs of cells)
% [faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] = getFaceAdjFromJnAdjGraph(edgeListInds,nodeEdges,...
%     junctionTypeListInds,jAnglesAll_alpha,boundaryEdges,edges2nodes,ws,edges2pixels);
disp('calculating adjacency graph of regions ...')
[faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] ...
    = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder,boundaryEdges);
% cellcogs = getCellCentroidsAll(setOfCells,edges2pixels,edgeListInds,...
%     sizeR,sizeC,edges2nodes,nodeInds);

[~,edgeOrientationsInds] = getEdgePriors(orientedScoreSpace3D,edges2pixels);
edgeOrientations = (edgeOrientationsInds-1).*orientationsStepSize;

% 
% cellStatesAll = getAllCellStates(setOfCells,edgeListInds,edgeOrientations,...
%     sizeR,sizeC,edges2pixels,nodeInds,edges2nodes);

% normalize input image
normalizedInputImage = imgIn./(max(max(imgIn)));

% cellPriors_old = getCellPriors_intensity(normalizedInputImage,setOfCells,edges2pixels,...
%     nodeInds,edges2nodes,cCell);
% get cell priors from RFC probability map
if ~exist('forest.mat','file')
    disp('RF for membrane classification not found. Training new classifier...')
    forest = trainRandomForest_pixelProb();
else
    load forest.mat
    disp('loaded pre-trained RF for membrane vs cell-interior classification')
end
regionPriors = regionScoreCalculator(forest,normalizedInputImage,setOfRegions,edges2pixels,...
    nodeInds,edges2nodes,cCell,wsIDsForRegions,ws);
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
            
% remove the edges connected to the boundaryNodes

% remove boundaryNodeEdgeListIDs from the offEdgeListIDs
offEdgeListIDs = setdiff(offEdgeListIDs,boundaryNodeEdgeListIDs);
          
% offEdgeListIDs = [];
% setting edgePriors
edgePriors(offEdgeListIDs) = offEdgeReward;
% visualize off edges
imgOffEdges = visualizeOffEdges(offEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC);
figure;imshow(imgOffEdges); title('visualization of edges turned off')
%% Backbone

onEdgeListIDs = getBackboneEdgeIDs(edgepixels,edgePriors,...
                lenThreshBB,priorThreshFracBB);
% removing boundary edges from the onEdgeListID list
onEdgeListIDs = setdiff(onEdgeListIDs,boundaryEdgeListInds);

bbNodeListInds = getJunctionsForEdges(edges2nodes,onEdgeListIDs);
bbNodePixInds = nodeInds(bbNodeListInds);
bbnodesVis = wsRegionBoundariesFromGraph;
bbnodesVis(bbNodePixInds) = 0.2;
figure;imagesc(bbnodesVis);

% visualize BB edges
imgBBEdges = visualizeOffEdges(onEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC);
figure;imshow(imgBBEdges); title('visualization of backbone 1')

edgePriors(onEdgeListIDs) = bbEdgeReward;
% test - enforcing some edges to be picked
clear onEdgeListIDs
onEdgeListIDs = [];

% visualize BB edges
imgBBEdges = visualizeOffEdges(onEdgeListIDs,edgepixels,nodeInds,sizeR,sizeC);
figure;imshow(imgBBEdges); title('visualization of backbone 2')
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
% [Aeq,beq] = getEqConstraints2(numEdges,jEdges,edges2pixels);
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
    model.vtype = 'B';  % binary variables
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
% get active edges and active nodes from x
ilpSegmentation = zeros(sizeR,sizeC);
% active edges
% consider the edgeID given in the first col of edges2pixels?? no need for
% this since we use edgepixels array which is already sans the skipped

% edges
onStateEdgeXind = 2:2:(numEdges*2);
onEdgeStates = x(onStateEdgeXind);
onEdgeInd = find(onEdgeStates>0.5);
offEdgeListInd = find(onEdgeStates<0.5);
onEdgePixelInds = getPixSetFromEdgeIDset(onEdgeInd,edgepixels);
offEdgePixelInds = getPixSetFromEdgeIDset(offEdgeListInd,edgepixels);
ilpSegmentation(onEdgePixelInds) = 1;

% active nodes 
fIndStop = 2*numEdges;
nodeInactiveStates_x = [];
nodeActivationVector = zeros(numel(nodeInds),1);    % stores 1 for active node list inds
nodeIndsActive = [];
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
        nodeInactiveStates_x = [nodeInactiveStates_x; inactiveness_nodes_i];
        activeStateNodeListInd = find(inactiveness_nodes_i<0.5);
        if(~isempty(activeStateNodeListInd))
            nodeListInd_i = junctionNodesListInds_i(activeStateNodeListInd);
            nodeActivationVector(nodeListInd_i) = 1;
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
            nodeIndsActive = [nodeIndsActive; nodeIndsActive_i];
        end
    end
end

inactiveNodePixInds = setdiff(nodeInds,nodeIndsActive);
offNodeIndList = find(ismember(nodeInds,inactiveNodePixInds));

activeContourPixels = find(ilpSegmentation);

% figure;imagesc(ilpSegmentation);title('ILP contours');
% reconstruct the edges with the values from the orientation filters (HSV)
% [output rgbimg] = reconstructHSVgauss_mv(orientedScoreSpace3D,orientations,...
%             barLength,barWidth,threshFrac,medianFilterH);
% get the active pixels
%output(:,:,3) = ilpSegmentation;
totX = numel(x);
numRegions = numel(regionPriors);
% get active foreground cells
regionStartPos = totX - numRegionVars + 2;
regionActivationVector = x(regionStartPos:2:totX);
activeRegionInd = find(regionActivationVector>0);
% get internal pixels for foreground cells
foregroundPixels = [];

% wsIDs = getWsIDsForCellIDs(ws,setOfCells,edges2pixels,nodeInds,...
%             edges2nodes,edgeListInds);

for i=1:numel(activeRegionInd)
    wsID_i = wsIDsForRegions(activeRegionInd(i));
    if(wsID_i==0)
        disp('problem with wsid check')
    else
%         cellPixInds_i = find(ws==wsID_i);
        cellPixInds_i = getInternalPixForCell(ws,wsID_i);
        foregroundPixels = [foregroundPixels; cellPixInds_i];
    end
end

% foregroundPixels = setdiff(foregroundPixels,offEdgePixelInds);

% make inactive edges inside foreground regions show as foreground
[inEdgeListInds,inEdgeIDs] = getInEdges(twoRegionEdges,regionActivationVector,...
                onEdgeStates,edges2regions,edgeListInds);
            
if(~isempty(inEdgeListInds))
    inEdgePixels = [];
    for i=1:numel(inEdgeListInds)
        inEdgePixels_i = edgepixels(inEdgeListInds(i),:);
        inEdgePixels_i = inEdgePixels_i(inEdgePixels_i>0);
        inEdgePixels = [inEdgePixels; inEdgePixels_i'];
    end
    foregroundPixels = [foregroundPixels; inEdgePixels];
end

% make inactive nodes inside foreground regions show as foreground
% nodeActivationVector = ~nodeInactiveStates_x;

inNodePixels = getInNodePixels(inEdgeIDs,nodeEdges,...
        nodeActivationVector,connectedJunctionIDs);
foregroundPixels = [foregroundPixels; inNodePixels];


foregroundPixels = setdiff(foregroundPixels,activeContourPixels);


% first visualize the contours
output_h = output(:,:,1);
output_s = output(:,:,2);
output_v = ilpSegmentation;

% create HSV image
% hsvImage = cat(3,output(:,:,1),output(:,:,2),ilpSegmentation);
hsvImage = cat(3,output_h,output_s,output_v);
% convert it to an RGB image
RGBimg = hsv2rgb(hsvImage);
% titleStr = sprintf('C = %d : lambda = %d',cNode,decayRate);
% titleStr = sprintf('C = %d',maxCost_direction);
figure;imshow(RGBimg)


% assign white to active (foreground) cells

output_h(foregroundPixels) = 1; 
output_s(foregroundPixels) = 0;
output_v(foregroundPixels) = 1;

% create HSV image
% hsvImage = cat(3,output(:,:,1),output(:,:,2),ilpSegmentation);
hsvImage_foreground = cat(3,output_h,output_s,output_v);
% convert it to an RGB image
RGBimg_foreground = hsv2rgb(hsvImage_foreground);
% titleStr = sprintf('C = %d : lambda = %d',cNode,decayRate);
% titleStr = sprintf('C = %d',maxCost_direction);
figure;imshow(RGBimg_foreground)
%% store extracted geometry in datastructures
[c_cellBorderEdgeIDs,c_cellBorderNodeIDs] = getCellBorderComponents(onEdgeInd,...
        edges2nodes,nodeEdges,edgeListInds);

cellBorderPixels = getCellBorderPixels(c_cellBorderEdgeIDs,...
            c_cellBorderNodeIDs,edgepixels,nodeInds,connectedJunctionIDs);
        
visualizeCellBorders = zeros(sizeR,sizeC);
visualizeCellBorders(cellBorderPixels) = 1;
figure;imshow(visualizeCellBorders)

% regions aggregating to form cells
% offEdgeIDList = edgeListInds(ismember(edgeListInds,offEdgeListInd)); 
offEdgeIDList = edgeListInds(offEdgeListInd);

[c_cells2regions,c_cellInternalEdgeIDs,c_cellIntNodeListInds] = getRegionsForOnCells(...
            faceAdj,activeRegionInd,offEdgeIDList,setOfRegions,wsIDsForRegions,ws,...
            offNodeIndList,edges2nodes,edgeListInds);
                
% visualize each cell in different colors
visualizeCells = zeros(sizeR,sizeC,3);
numCs = numel(c_cells2regions);
rMat = zeros(sizeR,sizeC);
gMat = zeros(sizeR,sizeC);
bMat = zeros(sizeR,sizeC);
for i=1:numCs
    % pick random color (RGB vals)
    R = rand(1); G = rand(1); B = rand(1);
    
    % get regions for this cell and the internal pixels. set RGB
    cellRegionList_i = c_cells2regions{i};
    regionPixels = getRegionPixels(cellRegionList_i,wsIDsForRegions,ws);
    rMat(regionPixels) = R;
    gMat(regionPixels) = G;
    bMat(regionPixels) = B;
    
    
    % get internal edges for this cell. set RGB.
    regionIntEdgeIDs = c_cellInternalEdgeIDs{i};
    regionIntEdgeListInds_logicalInd = ismember(edgeListInds,regionIntEdgeIDs);
    regionIntEdgePixels = edgepixels(regionIntEdgeListInds_logicalInd,:);
    regionIntEdgePixels = regionIntEdgePixels(regionIntEdgePixels>0);
    
    rMat(regionIntEdgePixels) = R;
    gMat(regionIntEdgePixels) = G;
    bMat(regionIntEdgePixels) = B;
      
    % get internal nodes for this cell. set RGB
    regionIntNodeListInds = c_cellIntNodeListInds{i};
    regionIntNodePixInds = getNodePixelsFromNodeInd...
                (regionIntNodeListInds,nodeInds,connectedJunctionIDs);
    
    rMat(regionIntNodePixInds) = R;
    gMat(regionIntNodePixInds) = G;
    bMat(regionIntNodePixInds) = B;
    
    visualizeCells(:,:,1) = rMat;
    visualizeCells(:,:,2) = gMat;
    visualizeCells(:,:,3) = bMat;
            
    
end

figure;imshow(visualizeCells);

segmentationOut = removeThickBorder(visualizeCells,marginSize);

