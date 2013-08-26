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

for i=1:length(trainingLabelImgNames)
  name = trainingLabelImgNames(i).name;
  im = imread(name);
  % get the edges that are inside the positive region
  edgepixels_i = cells_edgepixels_allTrainingImgs{i};
  numEdges_i = size(edgepixels_i,1);
  posPos = find(im(:,:,2)==255 & im(:,:,1)==0);
  posNeg = find(im(:,:,1)==255 & im(:,:,2)==0);
  for j=1:numEdges_i
    edgePix_edge_i = edgepixels(j,:);
    clear negativePix_j
    negativePix_j = posNeg(edgePix_edge_i);
    if(isempty(negativePix_j))
        % the edges is positive
        
    end
  end

  
  if ~isempty(posPos) || ~isempty(posNeg)
      load(strcat(name(1:LEN_IMG_IND),'_fm.mat'));
      fm = reshape(fm,size(fm,1)*size(fm,2),size(fm,3));
      fm(isnan(fm))=0;
      fmPos = [fmPos; fm(posPos,:)];
      fmNeg = [fmNeg; fm(posNeg,:)];
      clear fm;
  end
end

%% Train RFC

%% Visualization and evaluation of the trained RFC