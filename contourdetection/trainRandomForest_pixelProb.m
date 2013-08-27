function forest_pixelProb = trainRandomForest_pixelProb()

% Script to train Random Forest for pixel classification
%   0 - cell interior (bright) -> RED
%   1 - membrane (dark) -> GREEN

% Annotations for training images should be done only in green (+) and red
% (-)

%% parameters
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
% RF testing param

% param for membrane feature learning
  cs = 29;
  ms = 3;
  csHist = cs;

%% read training images
imgFiles_training = dir(pathForImages_training); % images for training

%% extract features. save separate feature matrices for each training image
% name of feature matrix: <imageName(1:ind)>_fm.mat
for i=1:length(imgFiles_training)
    name = imgFiles_training(i).name;
    % extract features only if feature matrix is not already presaved
    if ~exist(strcat(name(1:LEN_IMG_IND),'_fm.mat'),'file')
        disp('extracting features for the training images ...');
        im = norm01((imresize(imread(name),1)));
        % fm = getFeatures(im); % TODO
        fm  = membraneFeatures(im, cs, ms, csHist);
        save(strcat(name(1:LEN_IMG_IND),'_fm.mat'),'fm');
        clear fm
        clear im       
        disp('Feature extraction of training images done!')
        disp('Feature matrices saved as *_fm.mat files');
    else
        txt1 = sprintf('Found precalculated feature matrix for training image %d',i);
        disp(txt1);
    end
end

%% read training labels
trainingLabelImgNames = dir(pathForLabels_training);
fmPos = [];
fmNeg = [];

for i=1:length(trainingLabelImgNames)
  name = trainingLabelImgNames(i).name;
  im = imread(name);
  posPos = find(im(:,:,2)==255 & im(:,:,1)==0);
  posNeg = find(im(:,:,1)==255 & im(:,:,2)==0);
  
  if ~isempty(posPos) || ~isempty(posNeg)
      load(strcat(name(1:LEN_IMG_IND),'_fm.mat'));
      fm = reshape(fm,size(fm,1)*size(fm,2),size(fm,3));
      fm(isnan(fm))=0;
      fmPos = [fmPos; fm(posPos,:)];
      fmNeg = [fmNeg; fm(posNeg,:)];
      clear fm;
  end
end

clear posPos
clear posNeg

%% Training RFC
disp('Training ...')
disp('Number of training samples per class: ');
disp('dark (1): ');
disp(size(fmPos,1));
disp('bright (cell interior - 0): ');
disp(size(fmNeg,1));

y = [zeros(size(fmNeg,1),1);ones(size(fmPos,1),1)];
x = double([fmNeg;fmPos]);

extra_options.sampsize = [maxNumberOfSamplesPerClass, maxNumberOfSamplesPerClass];
if ~exist('forest_pixelProb.mat','file')
    forest_pixelProb = classRF_train(x, y, NUM_TREES,MTRY,extra_options);
    save forest_pixelProb.mat forest_pixelProb
else
    disp('forest_pixelProb.mat already exists!')
    load forest_pixelProb.mat
end

%% Evaluation of the learned forest
% Read test images
imgFiles_testing = dir(pathForImages_testing);
% calculate feature matrices for test images
for i=1:length(imgFiles_testing)
    name = imgFiles_testing(i).name;
    % extract features only if feature matrix is not already presaved
    if ~exist(strcat(name(1:LEN_IMG_IND),'_fm.mat'),'file')
        disp('extracting features for the test images ...');
        im = norm01((imresize(imread(name),1)));
        % fm = getFeatures(im); % TODO
        fm  = membraneFeatures(im, cs, ms, csHist);
        save(strcat(name(1:LEN_IMG_IND),'_fm.mat'),'fm');
        clear fm
        clear im       
        disp('Feature extraction of testing images done!')
        disp('Feature matrices saved as *_fm.mat files');
    end
end
% Read labels for test images
testingLabelImgNames = dir(pathForLabels_testing);

for i=1:length(testingLabelImgNames)
  name = testingLabelImgNames(i).name;
  im = imread(name);
  load(strcat(name(1:LEN_IMG_IND),'_fm.mat'));
  fm = reshape(fm,size(fm,1)*size(fm,2),size(fm,3));
  fm(isnan(fm))=0;
  clear fmNeg
  clear fmPos
  im=uint8Img(im(:,:,1));
  imsize = size(im);
  clear y
  clear im
  
  votes = zeros(imsize(1)*imsize(2),1);
  
  txt1 = sprintf('prediction for test image %d',i);
  disp(txt1);
  [y_h,v] = classRF_predict(double(fm), forest_pixelProb);
  txt1 = sprintf('prediction for test image %d done!',i);
  disp(txt1);
  votes = v(:,2);
  votes = reshape(votes,imsize);
  votes = double(votes)/max(votes(:));
  
  disp('visualization')

  im = imread(name);			% 
  % this illustration uses the thickened skeleton of the segmentation
  % this is the skeletonized view
  figure; imshow(makeColorOverlay(votes,im));
  imwrite(makeColorOverlay(votes,im),strcat(name(1:LEN_IMG_IND),'_overlay.tif'),'tif');
  % TODO: calculate pixel error
  
end
