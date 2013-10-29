function forest_edgeProb = trainRandomForest_edgeScore()

% OLDER version

% Script to train RFC for edge classification
%   0 - edge to be turned off (not part of contour)
%   1 - edge to be turned on (part of contour)

% Annotations for training images should be done only in green (1) and red
% (0)

%% Parameters
% use wild cards to allow for indices
pathForImages_training = '/home/thanuja/Dropbox/data/RF_training_edge/*_trainingImage.tif'; 
pathForImages_testing = '/home/thanuja/Dropbox/data/RF_training_edge/*_testingImg.tif';
pathForLabels_training = '/home/thanuja/Dropbox/data/RF_training_edge/*_trainingLabels.tif';
pathForLabels_testing = '/home/thanuja/Dropbox/data/RF_training_edge/*_testingLabels.tif';

maxNumberOfSamplesPerClass = 5000;
LEN_IMG_IND = 3;
% RF training param
NUM_TREES = 500;
MTRY = 5;

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

%% Read training images
imgFiles_training = dir(pathForImages_training); % images for training

%% Extract features for training images

fm = [];
numTrainingImgs = length(imgFiles_training); 
cells_edgepixels_allTrainingImgs = cell(1,numTrainingImgs);
disp('Edge feature extraction ... ');
for i=1:numTrainingImgs
    name = imgFiles_training(i).name;
    str1 = sprintf('image %d:',i);
    disp(str1)
    if ~exist(strcat(name(1:LEN_IMG_IND),'_efm.mat'),'file')
        imgIn0 = double(imread(name));
        % add border
        marginPixVal = min(min(imgIn0));
        imgIn = addThickBorder(imgIn0,marginSize,marginPixVal);
        [output,rgbimg,orientedScoreSpace3D] = getOFR(imgIn,orientations,...
                                barLength,barWidth,invertImg,threshFrac);
        OFR_mag = output(:,:,3);
        ws = watershed(OFR_mag);
        [sizeR,sizeC] = size(ws);
        disp('creating graph from watershed boundaries...');
        [adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs]...
            = getGraphFromWS(ws,output,0);
        edgepixels_i = edges2pixels;
        edgepixels_i(:,1) = []; % delete the first column which has edgeIDs
        cells_edgepixels_allTrainingImgs{i} = edgepixels_i;
        
        % extract features
        edgePrior = getEdgeUnaryAbs(edgepixels_i,OFR_mag);
        str2 = sprintf('calculating edge features for image %d',i);
        disp(str2)
        fm = edgeFeatures(imgIn,edgepixels_i,orientedScoreSpace3D,edgePrior);
        % save feature matrix
        save(strcat(name(1:LEN_IMG_IND),'_efm.mat'),'fm');
    else
        txt1 = sprintf('Found precalculated feature matrix for training image %d',i);
        disp(txt1);
    end
        
end

%% Read training labels
% load training labels
% 1 - green - membrane/mito
% 0 - red - cell interior

% for each edge check if it's in the on region. if yes, that is an on edge
trainingLabelImgNames = dir(pathForLabels_training);
fmPos = [];
fmNeg = [];

for i=1:numTrainingImgs
    
  name = trainingLabelImgNames(i).name;
  im = imread(name);
  % add border - TODO: what value to assign to the label border?
  marginPixVal = 0;
  im = addThickBorder(im,marginSize,marginPixVal);
  
  % get the edges that are inside the positive region
  edgepixels_i = cells_edgepixels_allTrainingImgs{i};
  numEdges_i = size(edgepixels_i,1);
  posPos = find(im(:,:,2)==255 & im(:,:,1)==0);
  posNeg = find(im(:,:,1)==255 & im(:,:,2)==0);
  posEdgeID_i = [];
  negEdgeID_i = [];
  for j=1:numEdges_i
    edgePix_edge_i = edgepixels_i(j,:);
    clear negativePix_j
    negativePix_j = posNeg(edgePix_edge_i);
    if(isempty(negativePix_j))
        % the edge is positive
        posEdgeID_i = [posEdgeID_i; j];
    else
        negEdgeID_i = [negEdgeID_i; j];
    end
  end

  if ~isempty(posPos) || ~isempty(posNeg)
      load(strcat(name(1:LEN_IMG_IND),'_efm.mat'));
      fm = reshape(fm,size(fm,1)*size(fm,2),size(fm,3));
      fm(isnan(fm))=0;
      fmPos = [fmPos; fm(posEdgeID_i,:)];
      fmNeg = [fmNeg; fm(negEdgeID_i,:)];
      clear fm;
  end
% visualize training data - on edges    
edgeMap = zeros(size(im));
positiveEdges = edgepixels_i(posEdgeID_i,:);
positiveEdgePixels = positiveEdges(positiveEdges>0);
edgeMap(positiveEdgePixels) = 1;
figure;imshow(edgeMap);title('training img - positive edges')
end

clear posPos
clear posNeg

%% Train RFC
disp('Training ...')
disp('Number of training samples per class: ');
disp('on edges: ');
disp(size(fmPos,1));
disp('off edges: ');
disp(size(fmNeg,1));

y = [zeros(size(fmNeg,1),1);ones(size(fmPos,1),1)];
x = double([fmNeg;fmPos]);

extra_options.sampsize = [maxNumberOfSamplesPerClass, maxNumberOfSamplesPerClass];
if ~exist('forest_edgeProb.mat','file')
    forest_edgeProb = classRF_train(x, y, NUM_TREES,MTRY,extra_options);
    disp('RFC learned for edge classification!')
    save forest_edgeProb.mat forest_edgeProb
    disp('saved forest forest_edgeProb.mat')
else
    disp('forest_edgeProb.mat already exists!')
    load forest_edgeProb.mat
end
%% Visualization and evaluation of the trained RFC
% Read test images
disp('Testing the learned RFC for edge classification ...')
imgFiles_testing = dir(pathForImages_testing);
% calculate feature matrices for test images
disp('Edge feature extraction for test images ...')
for i=1:length(imgFiles_testing)
    name = imgFiles_testing(i).name;
    % extract features only if feature matrix is not already presaved
    if ~exist(strcat(name(1:LEN_IMG_IND),'_efm.mat'),'file')
        imgIn0 = double(imread(name));
        % add border
        marginPixVal = min(min(imgIn0));
        imgIn = addThickBorder(imgIn0,marginSize,marginPixVal);
        [output,rgbimg,orientedScoreSpace3D] = getOFR(imgIn,orientations,...
                                barLength,barWidth,invertImg,threshFrac);
        OFR_mag = output(:,:,3);
        ws = watershed(OFR_mag);
        
        disp('creating graph from watershed boundaries...');
        [adjacencyMat,nodeEdges,edges2nodes,edges2pixels,connectedJunctionIDs] = getGraphFromWS(ws,output);
        edgepixels_i = edges2pixels;
        edgepixels_i(:,1) = []; % delete the first column which has edgeIDs

        % extract features
        edgePrior = getEdgeUnaryAbs(edgepixels_i,OFR_mag);
        str2 = sprintf('Calculating edge features for test image %d ...',i);
        disp(str2)
        fm = edgeFeatures(edgepixels_i,orientedScoreSpace3D,edgePrior);
        disp('done!')
        % save feature matrix
        save(strcat(name(1:LEN_IMG_IND),'_efm.mat'),'fm');
        % save edgepixels
        save(strcat(name(1:LEN_IMG_IND),'_edgepixels.mat'),'edgepixels_i');
    else
        txt1 = sprintf('Found precalculated feature matrix for test image %d',i);
        disp(txt1);
    end
    
end
% Read labels for test images
testingLabelImgNames = dir(pathForLabels_testing);

for i=1:length(testingLabelImgNames)
  name = testingLabelImgNames(i).name;
  im = imread(name);
  % add border
  marginPixVal = min(min(im));
  im = addThickBorder(im,marginSize,marginPixVal);
  [sizeR,sizeC] = size(im);
  load(strcat(name(1:LEN_IMG_IND),'_efm.mat'));
  load(strcat(name(1:LEN_IMG_IND),'_edgepixels.mat'));
  fm = reshape(fm,size(fm,1)*size(fm,2),size(fm,3));
  fm(isnan(fm))=0;
  clear fmNeg
  clear fmPos
  im=uint8Img(im(:,:,1));
  imsize = size(im);
  clear y
  clear im
  
  votes = zeros(imsize(1)*imsize(2),1);
  
  txt1 = sprintf('Classifying edges for test image %d ...',i);
  disp(txt1);
  [y_h,v] = classRF_predict(double(fm), forest_pixelProb);
  txt1 = sprintf('classifying edges for test image %d done!',i);
  disp(txt1);
  votes = v(:,2);
  votes = double(votes)/max(votes(:)); % on probability for each edge
  numEdges = size(votes,1);
  
  disp('visualization')
  edgeMap = zeros(sizeR,sizeC);
  for i=1:numEdges
      edgeMap(edgepixels_i(i)) = votes(i);
  end
  im = imread(name);			% 
  % this illustration uses the thickened skeleton of the segmentation
  % this is the skeletonized view
  figure; imshow(makeColorOverlay(edgeMap,im));
  imwrite(makeColorOverlay(edgeMap,im),strcat(name(1:LEN_IMG_IND),'_overlay.tif'),'tif');
  % TODO: calculate edge prediction error
  
end



