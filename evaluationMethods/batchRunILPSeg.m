% batch run of doILPseg

inputDir = '/home/thanuja/Dropbox/data/evaldata/input/';
outputDir = '/home/thanuja/Dropbox/data/evaldata/output/exp3/';
LEN_IMG_IND = 3;

% read files from directory
inputImageFileNames = strcat(inputDir,'/*_raw05.tif');
imgFiles_training = dir(inputImageFileNames); 
% doILP
for i=1:length(imgFiles_training)
    i
    name = imgFiles_training(i).name;
    % ILP
    inputImgFileName = strcat(inputDir,name);
    seg = doILP_w(inputImgFileName);
    % save
    saveFileName = strcat(name(1:LEN_IMG_IND),'_neuronseg05.png');
    saveFileName = strcat(outputDir,saveFileName);
    % save(saveFileName,'seg');    
    imwrite(seg,saveFileName,'png');
end

%% eval
trueLabelDir = '/home/thanuja/Dropbox/data/evaldata/labels/';
[VI, ME, SE, totNumNeur] = evaluateSegmentation(trueLabelDir,outputDir);