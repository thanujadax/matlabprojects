function trainRFCs3DILP()

createTrainingData = 1;
if(createTrainingData)
    disp('Training data creation start time:')
    datestr(clock, 0)
    createXYforRFtraining();
    disp('Training data creation stop time:')
    datestr(clock, 0)
end

% Train 3 random forest classifiers
disp('Training RFC start time:')
datestr(clock, 0)
%% File paths
% savePathRFC = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/RFCs';
% mkdir(savePathRFC);
% 
% fm_gridCells = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/fm/cellInteriorFMs/fm_cellInterior_sansBorderCells.mat';
% labels_gridCells = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/lblmat/gridCellLabels.mat';
% 
% fm_faces12 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/fm/cellFaceFMs/fm_faces12.mat';
% labels_faces12 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/lblmat/gridCellFace12Labels.mat';
% 
% fm_faces3456 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/fm/cellFaceFMs/fm_faces3456.mat';
% labels_faces3456 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/lblmat/gridCellFace3456Labels.mat';

savePathRFC = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/RFCs';
mkdir(savePathRFC);

fm_gridCells = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellInteriorFMs/fm_cellInterior_sansBorderCells.mat';
labels_gridCells = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellLabels.mat';

fm_faces12 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellFaceFMs/fm_faces12.mat';
labels_faces12 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellFace12Labels.mat';

fm_faces34 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellFaceFMs/fm_faces34.mat';
fm_faces56 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/cellFaceFMs/fm_faces56.mat';
labels_faces34 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellFace34Labels.mat';
labels_faces56 = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/lblmat/gridCellFace56Labels.mat';

forestName_1 = 'RFC_gridCells.mat';
forestName_2 = 'RFC_gridFaces12.mat';
forestName_3 = 'RFC_gridFaces3456.mat';
% forestName_3 = 'RFC_gridFaces34.mat'

%% RF training param
NUM_TREES = 500;
MTRY = 5;
maxNumberOfSamplesPerClass = 35000;
extra_options.sampsize = [maxNumberOfSamplesPerClass, maxNumberOfSamplesPerClass];
%% RFC1: gridCells
fm = importdata(fm_gridCells);
labels = importdata(labels_gridCells);
[x,y] = binaryClassBalancing(fm,labels);
disp('Training RFC for gridCell classification...')
RFC_gridCells = classRF_train(x, y, NUM_TREES,MTRY,extra_options);
disp('RFC trained for gridCell classification!')
saveFileName = fullfile(savePathRFC,forestName_1);
save(saveFileName,'RFC_gridCells','-v7.3');
disp('saved forest RFC_gridCells.mat')
%% RFC2: cell faces12 (xy)
fm = importdata(fm_faces12);
labels = importdata(labels_faces12);
[x,y] = binaryClassBalancing(fm,labels);
disp('Training RFC for gridFaces12 classification...')
RFC_gridFaces12 = classRF_train(x, y, NUM_TREES,MTRY,extra_options);
disp('RFC trained for gridFace12 classification!')
saveFileName = fullfile(savePathRFC,forestName_2);
save(saveFileName,'RFC_gridFaces12','-v7.3');
disp('saved forest RFC_gridFaces12.mat')

%% RFC3: cell faces3456 (xz and yz)
fm = importdata(fm_faces34);
fm2 = importdata(fm_faces56);
labels = importdata(labels_faces34);
labels2 = importdata(labels_faces56);
[x,y] = binaryClassBalancing([fm; fm2],[labels;labels2]);
disp('Training RFC for gridFaces3456 classification...')
RFC_gridFaces3456 = classRF_train(x, y, NUM_TREES,MTRY,extra_options);
disp('RFC trained for gridFace3456 classification!')
saveFileName = fullfile(savePathRFC,forestName_3);
save(saveFileName,'RFC_gridFaces3456','-v7.3');
disp('saved forest RFC_gridFaces3456.mat')

disp('Training RFC Finish time:')
datestr(clock, 0)