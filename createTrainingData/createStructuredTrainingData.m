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
[faceAdj,edges2regions,setOfRegions,twoRegionEdges,wsIDsForRegions] ...
    = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder,boundaryEdges);

%% get regions that matches individual neurons (connected components)

[labelImg_indexed,numLabels] = getLabelIndexImg(labelImage);
cells2regions = get

% get internal edges that bound two active regions
% get external edges that bound just one active region

