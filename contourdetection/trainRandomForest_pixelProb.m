function forest_pixelProb = trainRandomForest_pixelProb()

% Script to train Random Forest for pixel classification
%   0 - cell interior (bright) -> RED
%   1 - membrane (dark) -> GREEN

% Annotations for training images should be done only in green (+) and red
% (-)

%% parameters
pathForImages = ''; % use wild cards to allow for indices
pathForLabels = '';
maxNumberOfSamplesPerClass = 5000;
LEN_IMG_IND = 4;
% RF training params
NUM_TREES = 500;
MTRY = 5;

%% read images
imgFiles = dir(pathForImages);
labelFiles = dir(pathForLabels);

%% extract features. save separate feature matrices for each training image
% name of feature matrix: <imageName(1:6)>_fm.mat
for i=1:length(imgFiles)
    name = imgFiles(i).name;
    % extract features only if feature matrix is not already presaved
    if ~exist(strcat(name(1:LEN_IMG_IND),'_fm.mat'));
        disp('extracting membrane features');
        im = norm01((imresize(imread(imgNames(i).name),1)));
        fm = getFeatures(im);
        save(strcat(name(1:LEN_IMG_IND),'_fm.mat'),'fm');
        clear fm
        clear im       
    end
     
end

%% read training labels
trainingImgNames = dir(pathForLabels);
fmPos = [];
fmNeg = [];

for i=1:length(imgNames)
  name = imgNames(i).name;
  im = imread(name);
  posPos = find(im(:,:,2)==255 & im(:,:,1)==0);
  posNeg = find(im(:,:,1)==255 & im(:,:,2)==0);
  
  if length(posPos)>0 | length(posNeg)>0
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
forest_pixelProb = classRF_train(x, y, NUM_TREES,MTRY,extra_options);

