% evaluating trained RFCs
MAX_NUM_TEST_POINTS = 2000;

%% Inputs
% RFC_gridCells
% RFCfile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/RFCs/RFC_gridCells.mat';
% featureFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellInteriorFMs/fm_cellInterior_sansBorderCells.mat';
% labelFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellLabels.mat';

% RFC_faces12
% RFCfile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/RFCs/RFC_gridFaces12.mat';
% featureFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellFaceFMs/fm_faces12.mat';
% labelFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellFace12Labels.mat';

% RFC_faces34
RFCfile= '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/RFCs/RFC_gridFaces3456.mat';
featureFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellFaceFMs/fm_faces34.mat';
labelFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellFace34Labels.mat';

% RFC_faces56
% same RFC file as before
featureFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellFaceFMs/fm_faces56.mat';
labelFile = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellFace56Labels.mat';

%% Evaluation
% load RFC
RFC = importdata(RFCfile);
% get features
fm = importdata(featureFile);
% get labels
Y_tst = importdata(labelFile);
% predict
[Y_hat,v] = classRF_predict(fm,RFC);
% test error
err = length(find(Y_hat~=Y_tst))/length(Y_tst);

fprintf('\nexample 1: error rate %f\n', err);