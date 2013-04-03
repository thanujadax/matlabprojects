
% batch processing
savefilepath = '/home/thanuja/Dropbox/RESULTS/Jan/membraneOrientations_8_v2/';
% inputImagePath = '/home/thanuja/Dropbox/data/mitoData/07.tif';
inputImageDirPath = '/home/thanuja/Dropbox/data/membranes/';
% read tiff image
% inputImageStack = imread(inputImagePath,'png');
% numImg = size(inputImageStack,3);
numImg = 19; % just for testing
barLength = 15;
barWidth = 4;
invertImg = 0;

for i=0:numImg
    % for each image in the tiff image stack
    saveImgFileName = [num2str(i,'%02d') '.tif'];
    filename = [inputImageDirPath num2str(i,'%02d') '.tif'];
    % image_i = inputImageStack(:,:,i);
    image_i = imread(filename,'tif');
    [orientedScores,hsvOutput,rgbImg] = getOFR(image_i,barLength,barWidth,invertImg); 
    % save it in a folder
    % saveFileName = sprintf('HSVoutput%d.mat',i);
    % saveFileName = filenameOut;
    %saveCsvFileName = sprintf('HSVoutput%d.csv',i);
    %saveFileName = strcat(savefilepath,filenameOut);
    %saveCsvFileName = strcat(savefilepath,saveCsvFileName);
    %saveImgFileName = sprintf('HSVoutput%d.png',i);
    %save(saveFileName,'hsvOutput','orientedScores')
    %csvwrite(saveCsvFileName,hsvOutput);
    saveImgFileName = strcat(savefilepath,saveImgFileName);
%     titleStr = sprintf('rgbImage%d',i);
%     h = figure;imshow(rgbImg);title(titleStr);    
%     saveas(h,saveImgFileName);
    imwrite(rgbImg,saveImgFileName,'tif');
end