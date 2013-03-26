
% batch processing
savefilepath = '/home/thanuja/Dropbox/RESULTS/Jan/groundTruth1/';
inputImagePath = '/home/thanuja/Dropbox/data/mitoData/07.tif';
% inputImagePath = '/home/thanuja/Dropbox/data/mitoData/gettheorientations.png';
% read tiff image
inputImageStack = imread(inputImagePath,'tif');
% numImg = size(inputImageStack,3);
numImg = 1; % just for testing
barLength = 15;
barWidth = 4;
invertImg = 0;

for i=1:numImg
    % for each image in the tiff image stack
    image_i = inputImageStack(:,:,i);
    [orientedScores,hsvOutput,rgbImg] = getOFR(image_i,barLength,barWidth,invertImg); 
    % save it in a folder
    saveFileName = sprintf('HSVoutput%d.mat',i);
    %saveCsvFileName = sprintf('HSVoutput%d.csv',i);
    saveFileName = strcat(savefilepath,saveFileName);
    %saveCsvFileName = strcat(savefilepath,saveCsvFileName);
    saveImgFileName = sprintf('HSVoutput%d.png',i);
    save(saveFileName,'hsvOutput','orientedScores')
    %csvwrite(saveCsvFileName,hsvOutput);
    saveImgFileName = strcat(savefilepath,saveImgFileName);
%     titleStr = sprintf('rgbImage%d',i);
%     h = figure;imshow(rgbImg);title(titleStr);    
%     saveas(h,saveImgFileName);
    imwrite(rgbImg,saveImgFileName,'png');
end