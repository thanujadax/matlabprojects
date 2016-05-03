% function forestEdgeProb = trainRF_edgeProb()

% produces a random forest to classify edges of WS graph
%   1 - active edge (neuron boundary)
%   0 - inactive (neuron interior)

%% Parameters
% use wild cards to allow for indices

outputRoot = '/home/thanuja/projects/RESULTS/contours/20160503_edgeProbabilityRFC';

pathForImages_training = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/edgeProbability/raw/training'; 
pathForImages_testing = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/edgeProbability/raw/test';

pathForLabels_training = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/edgeProbability/labels/training';
pathForLabels_testing = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/edgeProbability/labels/test';

pathForMembranes_training = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/edgeProbability/membranes/training';
pathForMembranes_testing = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/edgeProbability/membranes/test';

fileNameString = '*.tif';

showIntermediate = 0;

maxNumberOfSamplesPerClass = 100000;
LEN_IMG_IND = 3;
% RF training param
NUM_TREES = 500;
MTRY = 7; % number of predictors sampled for splitting each tree

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
marginPixVal = 0;
%% read training data - x
rawImageFiles_training = dir(fullfile(pathForImages_training,fileNameString)); % raw images for training
labelImageFiles_training = dir(fullfile(pathForLabels_training,fileNameString)); % training labels
membraneFiles_training = dir(fullfile(pathForMembranes_training,fileNameString)); % membrane probabilities

numTrainingImgs = length(rawImageFiles_training);
% Extract edges 
% Extract features for edges
x = []; % feature matrix
y = []; % label matrix
disp('extracting edge features for training...')
for i=1:numTrainingImgs
    
    str1 = sprintf('training image %d:',i);
    disp(str1)
    
    rawImagePath_i = fullfile(pathForImages_training,rawImageFiles_training(i).name);
    labelImagePath_i = fullfile(pathForLabels_training,labelImageFiles_training(i).name);
    membraneProbMapPath_i = fullfile(pathForMembranes_training,membraneFiles_training(i).name);
    
    [c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,...
    c_extNodeInds,inactiveEdgeIDs,edgeListInds,edgepixels,OFR,edgePriors,OFR_mag,...
    boundaryEdgeIDs]...
        = createStructuredTrainingData(rawImagePath_i,labelImagePath_i);
    numEdgesTot = numel(edgeListInds);
    edgeListIndSequence = 1:numEdgesTot;
    % edges part of object boundaries - active
    activeEdgeIDs = getElementsFromCell(c_extEdgeIDs);
    % edges inside cells - inactive
    inactiveInternalEdgeIDs = getElementsFromCell(c_internalEdgeIDs);
    % append to the list of inactive edges outside the cells
    
    [numR,numC] = size(inactiveEdgeIDs);
    if(numC>1)
        inactiveEdgeIDs = inactiveEdgeIDs';
    end
    [numR,numC] = size(inactiveInternalEdgeIDs);
    if(numC>1)
        inactiveInternalEdgeIDs = inactiveInternalEdgeIDs';
    end
    
    inactiveEdgeIDs = [inactiveEdgeIDs; inactiveInternalEdgeIDs]; 
    
    % get their features
    [~,activeEdgeListInds_tr] = intersect(edgeListInds,activeEdgeIDs);
    [~,inactiveEdgeListInds_tr] = intersect(edgeListInds,inactiveEdgeIDs);
    edgeListInds_reordered_tr = [activeEdgeListInds_tr; inactiveEdgeListInds_tr];
    edgepixels_reordered_tr = edgepixels(edgeListInds_reordered_tr,:);
    edgePriors_reordered_tr = edgePriors(edgeListInds_reordered_tr);
    rawImage = double(imread(rawImagePath_i));
    rawImage = rawImage./(max(max(rawImage)));
    rawImage = addThickBorder(rawImage,marginSize,marginPixVal);
    
    membraneProbabilityMap = double(imread(membraneProbMapPath_i));
    membraneProbabilityMap = membraneProbabilityMap./(max(max(membraneProbabilityMap)));
    membraneProbabilityMap = addThickBorder(membraneProbabilityMap,marginSize,0);
    
    % clear fm
    fm = getEdgeFeatureMat(rawImage,edgepixels_reordered_tr,OFR,...
        edgePriors_reordered_tr,boundaryEdgeIDs,edgeListInds_reordered_tr,...
        membraneProbabilityMap);
    % append to the feature matrix x and the label matrix y
    x = [x; fm];
    numActiveEdges = numel(activeEdgeIDs);
    numInactiveEdges = numel(inactiveEdgeIDs);
    activeLabelVect = ones(numActiveEdges,1);
    inactiveLabelVect = zeros(numInactiveEdges,1);
    y = [y; activeLabelVect; inactiveLabelVect];
    
%     % debug code
%     checksum1 = size(fm,1) - numel(activeLabelVect) - numel(inactiveLabelVect);
%     if(checksum1~=0)
%         qq=1;
%     end
    
    
    % visualize active edges (red) and inactive edges (green)
    [sizeR,sizeC] = size(OFR_mag);
    visualizeR = zeros(sizeR,sizeC);
    visualizeG = zeros(sizeR,sizeC);
    visualizeB = zeros(sizeR,sizeC);

    activeEdgePixels_groundTruth = edgepixels(activeEdgeListInds_tr,:); 
    activeEdgePixels_groundTruth = activeEdgePixels_groundTruth(activeEdgePixels_groundTruth>0);
    visualizeR(activeEdgePixels_groundTruth) = 1;
    
    inactiveEdgePixels_groundTruth = edgepixels(inactiveEdgeListInds_tr,:); 
    inactiveEdgePixels_groundTruth = inactiveEdgePixels_groundTruth(inactiveEdgePixels_groundTruth>0);
    visualizeB(inactiveEdgePixels_groundTruth) = 1;
   
    visualization = zeros(sizeR,sizeC,3);
    visualization(:,:,1) = visualizeR;
    visualization(:,:,2) = visualizeG;
    visualization(:,:,3) = visualizeB;

    figure; imshow(visualization); title('training labels')
    
    
end
disp('feature extraction done! Saving feature matrix x and label vector y...')
% save x and y
save(fullfile(outputRoot,'x.mat'),'x');
save(fullfile(outputRoot,'y.mat'),'y');
disp('saved!')

%% Train RFC

totActiveEdges = sum(y==1);
totInactiveEdges = sum(y==0);
str1 = sprintf('Number of positive samples: %d', totActiveEdges);
disp(str1)
str1 = sprintf('Number of negative samples: %d', totInactiveEdges);
disp(str1)

disp('Training RF for edges...')
extra_options.sampsize = [maxNumberOfSamplesPerClass, maxNumberOfSamplesPerClass];
if ~exist('forestEdgeProb.mat','file')
    forestEdgeProb = classRF_train(x, y, NUM_TREES,MTRY,extra_options);
    disp('RFC learned for edge classification!')
    save(fullfile(outputRoot,'forestEdgeProbV7.mat'),forestEdgeProb,'-v7.3');
%     save forestEdgeProb.mat forestEdgeProb
    save(fullfile(outputRoot,'forestEdgeProb.mat'),forestEdgeProb);
    disp('saved forest forestEdgeProb.mat')
else
    disp('forestEdgeProb.mat already exists!')
    load forestEdgeProb.mat
end

clear x y

%% Testing, visualization and evaluation
% read test image
rawImageFiles_testing = dir(fullfile(pathForImages_testing,fileNameString));
labelImageFiles_testing = dir(fullfile(pathForLabels_testing,fileNameString)); % training labels
membraneFiles_testing = dir(fullfile(pathForMembranes_testing,fileNameString));
numTestingImgs = length(rawImageFiles_testing);

x0 = []; % feature matrix for test data
y0 = []; % treu label vector for test data

disp('extracting features for test images...')
activeEdgeListInds = [];
edgeListInds_reordered = [];

for i=1:numTestingImgs
    
    str1 = sprintf('test image %d:',i);
    disp(str1)
    
    rawImagePath_i = fullfile(pathForImages_testing,rawImageFiles_testing(i).name);
    labelImagePath_i = fullfile(pathForLabels_testing,labelImageFiles_testing(i).name);
    membraneProbMapPath_i = fullfile(pathForMembranes_testing,membraneFiles_testing(i).name);
    
    [c_cells2WSregions,c_internalEdgeIDs,c_extEdgeIDs,c_internalNodeInds,...
    c_extNodeInds,inactiveEdgeIDs,edgeListInds,edgepixels,OFR,edgePriors,OFR_mag,...
    boundaryEdgeIDs]...
        = createStructuredTrainingData(rawImagePath_i,labelImagePath_i);
    % edges part of object boundaries - active
    activeEdgeIDs = getElementsFromCell(c_extEdgeIDs);
    % edges inside cells - inactive
    inactiveInternalEdgeIDs = getElementsFromCell(c_internalEdgeIDs);
    % append to the list of inactive edges outside the cells
    inactiveEdgeIDs = [inactiveEdgeIDs; inactiveInternalEdgeIDs]; 
    inactiveEdgeIDs = unique(inactiveEdgeIDs);
    
    % get their features
    [~,activeEdgeListInds] = intersect(edgeListInds,activeEdgeIDs);
    [~,inactiveEdgeListInds] = intersect(edgeListInds,inactiveEdgeIDs);
    edgeListInds_reordered = [activeEdgeListInds; inactiveEdgeListInds];
    edgepixels_reordered = edgepixels(edgeListInds_reordered,:);
    edgePriors_reordered = edgePriors(edgeListInds_reordered);
    
    rawImage = double(imread(rawImagePath_i));
    rawImage = rawImage./(max(max(rawImage)));
    rawImage = addThickBorder(rawImage,marginSize,marginPixVal);
    
    membraneProbabilityMap = double(imread(membraneProbMapPath_i));
    membraneProbabilityMap = membraneProbabilityMap./(max(max(membraneProbabilityMap)));
    membraneProbabilityMap = addThickBorder(membraneProbabilityMap,marginSize,0);
    
    clear fm
    fm = getEdgeFeatureMat(rawImage,edgepixels_reordered,OFR,...
        edgePriors_reordered,boundaryEdgeIDs,edgeListInds_reordered,...
        membraneProbabilityMap);
    % append to the feature matrix x and the label matrix y
    x0 = [x0; fm];
    numActiveEdges = numel(activeEdgeIDs);
    numInactiveEdges = numel(inactiveEdgeIDs);
    activeLabelVect = ones(numActiveEdges,1);
    inactiveLabelVect = zeros(numInactiveEdges,1);
    y0 = [y0; activeLabelVect; inactiveLabelVect];
    
end

disp('done. saving feature matrix for training data')
save(fullfile(outputRoot,'x0.mat'),'x0');
save(fullfile(outputRoot,'y0.mat'),'y0');
% predicting labels for test data
[y_h,v] = classRF_predict(double(x0), forestEdgeProb);
predictedEdgeProbabilities = v(:,2);
predictedEdgeProbabilities = predictedEdgeProbabilities./(max(predictedEdgeProbabilities));
% visualize
% active edges - ground truth - green
% active edges - prediction - red
[sizeR,sizeC] = size(OFR_mag);
visualizeR = zeros(sizeR,sizeC);
visualizeG = zeros(sizeR,sizeC);
visualizeB = zeros(sizeR,sizeC);

activeEdgePixels_groundTruth = edgepixels(activeEdgeListInds,:); 
activeEdgePixels_groundTruth = activeEdgePixels_groundTruth(activeEdgePixels_groundTruth>0);
visualizeG(activeEdgePixels_groundTruth) = 1;

activeEdgeListInds_predicted = edgeListInds_reordered(logical(y_h));
activeEdgePixels_predicted = edgepixels(activeEdgeListInds_predicted,:);
activeEdgePixels_predicted = activeEdgePixels_predicted(activeEdgePixels_predicted>0);
visualizeR(activeEdgePixels_predicted) = 1;

visualization = zeros(sizeR,sizeC,3);
visualization(:,:,1) = visualizeR;
visualization(:,:,2) = visualizeG;
visualization(:,:,3) = visualizeB;

figure; imshow(visualization)

% prediction error
numTestEdges = numel(y_h)
predictionError = (sum(abs(y0 - y_h)))/numTestEdges

% visualize edge probabilities
visualizeEdgeProbPred = zeros(sizeR,sizeC);
visualizeEdgePriorProb = zeros(sizeR,sizeC);
numEdges = numel(edgeListInds_reordered);
for i=1:numEdges
    edgePixels_i = edgepixels(edgeListInds_reordered(i),:);
    edgePixels_i = edgePixels_i(edgePixels_i>0);
    visualizeEdgeProbPred(edgePixels_i) = predictedEdgeProbabilities(i);
    visualizeEdgePriorProb(edgePixels_i) = edgePriors_reordered(i);
end
figure;imagesc(visualizeEdgeProbPred);title('predicted probabilities')
figure;imagesc(visualizeEdgePriorProb);title('prior probabilities')

imwrite(visualizeEdgeProbPred,fullfile(outputRoot,'prediction.tif'),'tif');