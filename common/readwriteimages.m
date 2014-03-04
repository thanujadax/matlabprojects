% Read png file
% filename = '/home/thanuja/matlabprojects/data/mitoData/labels/stem1_membranes.png';
% filename = '/home/thanuja/matlabprojects/data/mitoData/stem.tiff';
filename = '/home/thanuja/Dropbox/data/em_2013january/neurons/07.tiff';
% filename = '/home/thanuja/Dropbox/data/RF_training_mem/I00_trainingLabels.tif';
% filename = '/home/thanuja/Dropbox/data/em_2013january/neurons/01.tiff';
fmt = 'tif';
A = imread(filename, fmt);
% C = rgb2gray(A);
% write to another file
% writefile = '/home/thanuja/Dropbox/data/em_2013january/samples/raw00_512.png';
% writefile = '/home/thanuja/Dropbox/data/RF_training_edge/I01_trainingLabels.tif';
writeFilePath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/neuron/';

dimx = 512;
dimy = 512;

startRow = 1;
stopRow = dimy;

startCol = 1;
stopCol = dimx;

numDim = 3;
writeFileName = '07.png';
writeFileName = strcat(writeFilePath,writeFileName);
B = A(startRow:stopRow,startCol:stopCol,:);
imwrite(B,writeFileName,'png')
figure;imshow(B);

% k = 00; % file index
% for i=1:4
%     for j=1:4
%         B = A(startRow:stopRow,startCol:stopCol,:);
%         writeName = sprintf('I%02d_trainingLabels.tif',k);
%         writeFileName = strcat(writeFilePath,writeName);
%         disp(writeFileName)
%         imwrite(B,writeFileName,'tif')
%         % figure;imshow(B);
%         
%         startCol = stopCol + 1;
%         stopCol = stopCol + dimx;
%         k = k + 1;
%     end
%     startRow = stopRow + 1;
%     stopRow = stopRow + dimy;
%     startCol = 1;
%     stopCol = dimx;
% end
        
