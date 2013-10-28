function forestEdgeProb = trainRF_edgeProb()

% produces a random forest to classify edges of WS graph
%   1 - active edge (neuron boundary)
%   0 - inactive (neuron interior)

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

%% read training data - x
imgFiles_training = dir(pathForImages_training); % images for training
% Extract edges
[internalEdgePixels, extEdgePixels] = createStructuredTrainingData(rawImagePath,labelImagePath);

% Extract features for edges


% read training labels - y



%% Train RFC


%% Visualization and evaluation



