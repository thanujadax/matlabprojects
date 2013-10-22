% evaluation script

function [VI,RI,ARI] = evaluateSegmentation()


trueLabelDir = '/home/thanuja/Dropbox/data/evaldata/labels/';
segDir = '/home/thanuja/Dropbox/data/evaldata/output/';

LEN_IMG_IND = 3;

% read files from directory
inputLabelFileNames = '/home/thanuja/Dropbox/data/evaldata/labels/*.tif';
imgFiles_labels = dir(inputLabelFileNames); 

inputSegFileNames = '/home/thanuja/Dropbox/data/evaldata/output/*.png';
imgFiles_seg = dir(inputSegFileNames); 

numImages = length(imgFiles_labels);

% ri = zeros(numImages,1);
% gce = zeros(numImages,1);
VI = zeros(numImages,1);
RI = zeros(numImages,1);
AR = zeros(numImages,1);


% evaluation metrics
for i=1:length(imgFiles_labels)
    i
    nameLabel = imgFiles_labels(i).name;
    nameSeg = imgFiles_seg(i).name;
    
    labelFileName = strcat(trueLabelDir,nameLabel);
    segFileName = strcat(segDir,nameSeg);
    
    imLabel = imread(labelFileName);
    imSeg = imread(segFileName);
    
    % convert into integer labels
    imLabel_int = im2uint16(rgb2gray(imLabel));
    imSeg_int = im2uint16(rgb2gray(imSeg));
    
    [L1m,L1n] = size(imLabel_int);
    [L2m,L2n] = size(imSeg_int);

    numL1 = L1m*L1n;
    numL2 = L2m*L2n;

    label_vec = reshape(imLabel_int,numL1,1);
    seg_vec = reshape(imSeg_int,numL2,1);
    
    % [ri(i),gce(i),vi(i)]=compare_segmentations(imLabel_int,imSeg_int);
    % vi(i) = varinfo(imLabel_int,imSeg_int);
    [variationOfInf,~] = vi(label_vec,seg_vec);
    VI(i) = variationOfInf;
    
    [AR(i),RI(i),~,~]=RandIndex(label_vec,seg_vec);
    
end

rimat = strcat(segDir,'ri.mat');
arimat = strcat(segDir,'ari.mat');
% gcemat = strcat(segDir,'gce.mat');
vimat = strcat(segDir,'vi.mat');

% save(rimat,'ri');
% save(gcemat,'gce');
save(vimat,'VI');
save(rimat,'RI');
save(arimat,'AR');
