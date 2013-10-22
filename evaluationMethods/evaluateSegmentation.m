% evaluation script

function [ri,gce,vi] = evaluateSegmentation()


trueLabelDir = '/home/thanuja/Dropbox/data/evaldata/labels/';
segDir = '/home/thanuja/Dropbox/data/evaldata/output/';

LEN_IMG_IND = 3;

% read files from directory
inputLabelFileNames = '/home/thanuja/Dropbox/data/evaldata/labels/*.tif';
imgFiles_labels = dir(inputLabelFileNames); 

inputSegFileNames = '/home/thanuja/Dropbox/data/evaldata/output/*.png';
imgFiles_seg = dir(inputSegFileNames); 

numImages = length(imgFiles_labels);

ri = zeros(numImages,1);
gce = zeros(numImages,1);
vi = zeros(numImages,1);

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
    
    [ri(i),gce(i),vi(i)]=compare_segmentations(imLabel_int,imSeg_int);
    rimat = strcat(segDir,'ri.mat');
    gcemat = strcat(segDir,'gce.mat');
    vimat = strcat(segDir,'vi.mat');
    
    save(rimat,'ri');
    save(gcemat,'gce');
    save(vimat,'vi');
end