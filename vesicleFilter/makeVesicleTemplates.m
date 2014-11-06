% script to extract vesicle templates
function makeVesicleTemplates()
imageFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/smallSet1/1.tif';
templatesPath = '/home/thanuja/Dropbox/data/fibsem/smallset1/vesicles/vesicleFilter2D/templates1';
templateFileName = 't1.png';
templateFullFile = fullfile(templatesPath,templateFileName);

% read image
inputImageMat = double(imread(imageFileName));
inputImageMat = inputImageMat ./255;
figure;imshow(inputImageMat);
% normalization (get rid of DC comp = rescale between -1 and +1)

% specify template corners
rstart = 0;
rstop = 10;
cstart = 0;
cstop = 10;

% extract template
outputTemplate = extractTemplate(inputImageMat,rstart,rstop,cstart,cstop);
figure; imagesc(outputTemplate)

% save template
imwrite(outputTemplate,templateFullFile,'png');
