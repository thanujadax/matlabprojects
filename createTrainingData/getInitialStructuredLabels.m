function [c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,...
    c_extNodeInds,edgeListInds,edgepixels,nodeInds,ws,...
            inactiveEdgeLIDs,offWsIDs,setOfRegions,edges2regions,edges2nodes,...
            nodeEdges,rawImage,labelImage,strDataVisualization,...
            edgeLID_confVector,wsID_confVector,...
            nodeLID_confVector] = ...
    getInitialStructuredLabels(rawImagePath,labelImagePath)
% create structured training labels: edges, nodes and regions
% version 2. 2013.12.05
% uses directed edges instead of undirected edges

% Inputs:
%   raw EM image
%   labels for the raw image: neurons, mitochondria? (later)

% Outputs:
%   for each cell, contains
%           extEdgeLIDs
%           internalEdgeIDs
%           extNodeLIDs
%           internalNodeLIDs
%           wsIDs

%% Parameters
showIntermediate = 0;

orientationsStepSize = 10;
orientations = 0:orientationsStepSize:350;

barLength = 13;     % should be odd
barWidth = 4;       %
b_imWithBorder = 1; % 1 to add thick border
marginSize = ceil(barLength/2);
marginPixVal = 0;
addBorder = ceil(barLength/2);
threshFrac = 0.1;
medianFilterH = 0;
invertImg = 1;      % 1 for EM images when input image is taken from imagePath

% rawImagePath = '/home/thanuja/Dropbox/data/evaldata/input/I03_raw05.tif';
% labelImagePath = '/home/thanuja/Dropbox/data/evaldata/labels2/I03_neuronLabels05.tif';

%% Read inputs. Perform initial processing
rawImage = double(imread(rawImagePath));
rawImage = rawImage./(max(max(rawImage)));
labelImage = imread(labelImagePath);

% add thick border
if(b_imWithBorder)
    rawImage = addThickBorder(rawImage,marginSize,marginPixVal);
    labelImage = addThickBorder(labelImage,marginSize,marginPixVal);
end
disp('Inputs read')
%% Extract WS graph for raw image
[HSVmat,rgbimg,OFR] = getOFR(rawImage,orientations,...
                            barLength,barWidth,invertImg,threshFrac);
OFR_mag = HSVmat(:,:,3);    % OFR value matrix

ws = watershed(OFR_mag);
[sizeR,sizeC] = size(ws);

disp('Calculating ws graph...')
[adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs,selfEdgePixelSet] ...
    = getGraphFromWS(ws,HSVmat,showIntermediate);
disp('done')
junctionTypeListInds = getJunctionTypeListInds(nodeEdges);
nodeInds = nodeEdges(:,1);                  % indices of the junction nodes
edgeListInds = edges2pixels(:,1);

clusterNodeIDs = connectedJunctionIDs(:,1); % indices of the clustered junction nodes

wsRegionBoundariesFromGraph = zeros(sizeR,sizeC);
wsRegionBoundariesFromGraph(nodeInds) = 0.7;          % junction nodes
wsRegionBoundariesFromGraph(clusterNodeIDs) = 0.5;    % cluster nodes
[nre,nce] = size(edges2pixels);  % first column is the edgeID
edgepixels = edges2pixels(:,2:nce);
wsRegionBoundariesFromGraph(edgepixels(edgepixels>0)) = 1; % edge pixels

boundaryEdgeIDs = getBoundaryEdges2(wsRegionBoundariesFromGraph,barLength,edgepixels,...
    nodeEdges,edgeListInds,showIntermediate);

disp('Calculating face adjacency graph...')
[faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] ...
    = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder,boundaryEdgeIDs);
disp('done')

% edgePriors = getEdgeUnaryAbs(edgepixels,OFR_mag);
%% get regions that matches individual neurons (connected components)

[labelImg_indexed,numLabels] = getLabelIndexImg(labelImage);

disp('Calculating initial labels using ground truth segmentations...')
[c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,c_extNodeInds]...
            = getCells2WSregions(labelImg_indexed,ws,numLabels,setOfRegions,...
            edgeListInds,edges2nodes);
disp('done')
% get offWsIDs
allCellWsIDs = getElementsFromCell(c_cells2WSregions);
numWsRegions = max(max(ws));
wsRegionSequence = 1:numWsRegions;
offWsIDs = setdiff(wsRegionSequence,allCellWsIDs);

numEdges = numel(edgeListInds);
edgeLIDsequence = 1:numEdges;

inactiveEdgeLIDs = edgeLIDsequence;     % initializing the list of inactive edgeIDs outside cells

%% Visualization

% visualize internal and external edges
strDataVisualization = zeros(sizeR,sizeC,3);
edgeMat2D_R = zeros(sizeR,sizeC);
edgeMat2D_G = zeros(sizeR,sizeC);
edgeMat2D_B = zeros(sizeR,sizeC);

for i=1:numLabels
    R_i = rand(1);
    G_i = rand(1);
    B_i = rand(1);
    
    R_e = rand(1);
    G_e = rand(1);
    B_e = rand(1);
    
    R_in = rand(1);
    G_in = rand(1);
    B_in = rand(1);
    
    R_en = rand(1);
    G_en = rand(1);
    B_en = rand(1);
    
    internalEdgeIDs_i = c_internalEdgeIDs{i};
    [~,internalEdgeLIDs_i] = intersect(edgeListInds,internalEdgeIDs_i);
    
    % update inactive edgeID list
    inactiveEdgeLIDs = setdiff(inactiveEdgeLIDs,internalEdgeLIDs_i);
    
    internalEdgeListInds_logical = ismember(edgeListInds,internalEdgeIDs_i);
    internalEdgePixels = edgepixels(internalEdgeListInds_logical,:);
    internalEdgePixels = internalEdgePixels(internalEdgePixels>0);
    
    edgeMat2D_R(internalEdgePixels) = R_i;
    edgeMat2D_G(internalEdgePixels) = G_i;
    edgeMat2D_B(internalEdgePixels) = B_i;
    
    extEdgeIDs_i = c_extEdgeIDs{i};
    [~,extEdgeLIDs_i] = intersect(edgeListInds,extEdgeIDs_i);
    
    % update inactive edgeID list
    inactiveEdgeLIDs = setdiff(inactiveEdgeLIDs,extEdgeLIDs_i);
    
    extEdgeListInds_logical = ismember(edgeListInds,extEdgeIDs_i);
    extEdgePixels = edgepixels(extEdgeListInds_logical,:);
    extEdgePixels = extEdgePixels(extEdgePixels>0);
    
    edgeMat2D_R(extEdgePixels) = R_e;
    edgeMat2D_G(extEdgePixels) = G_e;
    edgeMat2D_B(extEdgePixels) = B_e;
    
    internalNodeListInds = c_internalNodeInds{i};
    intNodePixInds = getNodePixelsFromNodeInd...
                (internalNodeListInds,nodeInds,connectedJunctionIDs);
            
    edgeMat2D_R(intNodePixInds) = R_in;
    edgeMat2D_G(intNodePixInds) = G_in;
    edgeMat2D_B(intNodePixInds) = B_in;
            
    extNodeListInds = c_extNodeInds{i};
    extNodePixInds = getNodePixelsFromNodeInd...
                (extNodeListInds,nodeInds,connectedJunctionIDs);
    
    edgeMat2D_R(extNodePixInds) = R_en;
    edgeMat2D_G(extNodePixInds) = G_en;
    edgeMat2D_B(extNodePixInds) = B_en;
            
end

% visualize inactive edges outside cells
% inactiveEdgeListInds_logical = ismember(edgeListInds,inactiveEdgeIDs);
inactiveEdgePixels = edgepixels(inactiveEdgeLIDs,:);
inactiveEdgePixels = inactiveEdgePixels(inactiveEdgePixels>0);
R_ina = rand(1);
G_ina = rand(1);
B_ina = rand(1);
edgeMat2D_R(inactiveEdgePixels) = R_ina;
edgeMat2D_G(inactiveEdgePixels) = G_ina;
edgeMat2D_B(inactiveEdgePixels) = B_ina;


% Visualize the entire thing
strDataVisualization(:,:,1) = edgeMat2D_R;
strDataVisualization(:,:,2) = edgeMat2D_G;
strDataVisualization(:,:,3) = edgeMat2D_B;

% figure;imshow(strDataVisualization)

