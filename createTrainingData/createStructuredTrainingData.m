% function trainingData = createStructuredTrainingData()
% create structured training labels: edges, nodes and regions

% Inputs
%   raw EM image
%   labels for the raw image: neurons, mitochondria? (later)

%% Parameters
showIntermediate = 1;

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

rawImagePath = '/home/thanuja/Dropbox/data/evaldata/input/I03_raw05.tif';
labelImagePath = '/home/thanuja/Dropbox/data/evaldata/labels2/I03_neuronLabels05.tif';

%% Read inputs. Perform initial processing
rawImage = imread(rawImagePath);
labelImage = imread(labelImagePath);

% add thick border
if(b_imWithBorder)
    rawImage = addThickBorder(rawImage,marginSize,marginPixVal);
    labelImage = addThickBorder(labelImage,marginSize,marginPixVal);
end
%% Extract WS graph for raw image
[HSVmat,rgbimg,orientedScoreSpace3D] = getOFR(rawImage,orientations,...
                            barLength,barWidth,invertImg,threshFrac);
OFR_mag = HSVmat(:,:,3);    % OFR value matrix

ws = watershed(OFR_mag);
[sizeR,sizeC] = size(ws);

[adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs,selfEdgePixelSet] ...
    = getGraphFromWS(ws,HSVmat,showIntermediate);

nodeInds = nodeEdges(:,1);                  % indices of the junction nodes
edgeListInds = edges2pixels(:,1);

clusterNodeIDs = connectedJunctionIDs(:,1); % indices of the clustered junction nodes

wsRegionBoundariesFromGraph = zeros(sizeR,sizeC);
wsRegionBoundariesFromGraph(nodeInds) = 0.7;          % junction nodes
wsRegionBoundariesFromGraph(clusterNodeIDs) = 0.5;    % cluster nodes
[nre,nce] = size(edges2pixels);  % first column is the edgeID
edgepixels = edges2pixels(:,2:nce);
wsRegionBoundariesFromGraph(edgepixels(edgepixels>0)) = 1; % edge pixels

boundaryEdges = getBoundaryEdges2(wsRegionBoundariesFromGraph,barLength,edgepixels,...
    nodeEdges,edgeListInds,showIntermediate);

[faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] ...
    = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder,boundaryEdges);

%% get regions that matches individual neurons (connected components)

[labelImg_indexed,numLabels] = getLabelIndexImg(labelImage);
[c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,c_extNodeInds]...
            = getCells2WSregions(labelImg_indexed,ws,numLabels,setOfRegions,...
            nodeEdges,connectedJunctionIDs);
% visualize internal and external edges
edgeVisualization = zeros(sizeR,sizeC,3);
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
    
    internalEdgeIDs_i = c_internalEdgeIDs{i};
    internalEdgeListInds_logical = ismember(edgeListInds,internalEdgeIDs_i);
    internalEdgePixels = edgepixels(internalEdgeListInds_logical,:);
    internalEdgePixels = internalEdgePixels(internalEdgePixels>0);
    
    edgeMat2D_R(internalEdgePixels) = R_i;
    edgeMat2D_G(internalEdgePixels) = G_i;
    edgeMat2D_B(internalEdgePixels) = B_i;
    
    extEdgeIDs_i = c_extEdgeIDs{i};
    extEdgeListInds_logical = ismember(edgeListInds,extEdgeIDs_i);
    extEdgePixels = edgepixels(extEdgeListInds_logical,:);
    extEdgePixels = extEdgePixels(extEdgePixels>0);
    
    edgeMat2D_R(extEdgePixels) = R_e;
    edgeMat2D_G(extEdgePixels) = G_e;
    edgeMat2D_B(extEdgePixels) = B_e;
        
end
edgeVisualization(:,:,1) = edgeMat2D_R;
edgeVisualization(:,:,2) = edgeMat2D_G;
edgeVisualization(:,:,3) = edgeMat2D_B;

figure;imshow(edgeVisualization)

