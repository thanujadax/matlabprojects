function forest_edgeScore = trainRandomForest_edgeScore()

% Script to train RFC for edge classification
%   0 - edge to be turned off (not part of contour)
%   1 - edge to be turned on (part of contour)

% Annotations for training images should be done only in green (1) and red
% (0)

%% Parameters
% use wild cards to allow for indices
pathForImages_training = '/home/thanuja/Dropbox/data/RF_training_mem/*_trainingImage.tif'; 
pathForImages_testing = '/home/thanuja/Dropbox/data/RF_training_mem/*_testingImg.tif';
pathForLabels_training = '/home/thanuja/Dropbox/data/RF_training_mem/*_trainingLabels.tif';
pathForLabels_testing = '/home/thanuja/Dropbox/data/RF_training_mem/*_testingLabels.tif';

maxNumberOfSamplesPerClass = 5000;
LEN_IMG_IND = 3;
% RF training param
NUM_TREES = 500;
MTRY = 5;

orientations = 0:10:350;
orientationsStepSize = 10;

barLength = 13; % should be odd
barWidth = 4; %
marginSize = ceil(barLength/2);
% marginPixVal = 0.1;
addBorder = ceil(barLength/2);
threshFrac = 0.1;
medianFilterH = 0;
invertImg = 1;      % 1 for EM images when input image is taken from imagePath

%% Read training images
imgFiles_training = dir(pathForImages_training); % images for training

%% Extract features for training images

fm = [];
numTrainingImgs = length(imgFiles_training); 
cells_edgepixels_allTrainingImgs = cell(1,numTrainingImgs);
for i=1:numTrainingImgs
    imgIn0 = double(imread(imgFiles_training(i).name));
    % add border
    marginPixVal = min(min(imgIn0));
    imgIn = addThickBorder(imgIn0,marginSize,marginPixVal);
    [output,rgbimg,orientedScoreSpace3D] = getOFR(imgIn,...
                            barLength,barWidth,invertImg,threshFrac);
    OFR_mag = output(:,:,3);
    ws = watershed(OFR_mag);
    [sizeR,sizeC] = size(ws);
    disp('creating graph from watershed boundaries...');
    [adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs] = getGraphFromWS(ws,output);
    edgepixels_i = edge2pixels;
    edgepixels_i(:,1) = []; % delete the first column which has edgeIDs
    cells_edgepixels_allTrainingImgs{i} =  edgepixels_i;
    % extract features
    edgePrior = getEdgeUnaryAbs(edgepixels,OFR_mag);
    fm_i = edgeFeatures(edgepixels_i,orientedScoreSpace3D,edgePrior);
    fm = [fm;fm_i];
end

%% Read training labels
% load training labels
% 1 - green - membrane/mito
% 0 - red - cell interior

% for each edge check if it's in the on region. if yes, that is an on edge
trainingLabelImgNames = dir(pathForLabels_training);
fmPos = [];
fmNeg = [];
cell_posEdgeIDs_all = cell(1,numTrainingImgs);
cell_negEdgeIDs_all = cell(1,numTrainingImgs);

totPosEdges = 0;
totNegEdges = 0;

for i=1:numTrainingImgs
  name = trainingLabelImgNames(i).name;
  im = imread(name);
  % get the edges that are inside the positive region
  edgepixels_i = cells_edgepixels_allTrainingImgs{i};
  numEdges_i = size(edgepixels_i,1);
  posPos = find(im(:,:,2)==255 & im(:,:,1)==0);
  posNeg = find(im(:,:,1)==255 & im(:,:,2)==0);
  posEdgeID_i = [];
  negEdgeID_i = [];
  for j=1:numEdges_i
    edgePix_edge_i = edgepixels(j,:);
    clear negativePix_j
    negativePix_j = posNeg(edgePix_edge_i);
    if(isempty(negativePix_j))
        % the edge is positive
        posEdgeID_i = [posEdgeID_i; j];
    else
        negEdgeID_i = [negEdgeID_i; j];
    end
  end
  
  totPosEdges = totPosEdges + numel(posEdgeID_i); 
  totNegEdges = totNegEdges + numel(negEdgeID_i);
  
  cell_posEdgeIDs_all{i} = posEdgeID_i;
  cell_negEdgeIDs_all{i} = negEdgeID_i;
  
  % TODO: build x and y within this loop

end

% create input training data matrix x
totNumEdges = totPosEdges + totNegEdges;
x = double


%% Train RFC

%% Visualization and evaluation of the trained RFC