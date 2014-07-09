% Read file

% filename = '/home/thanuja/projects/drosophila-l3/stack2/raw/00.tif';
% filename = '/home/thanuja/projects/drosophila-l3/stack2/classification/schmidhuber/median_filtered/neurons/neurons0000.png';
filename = '/home/thanuja/projects/drosophila-l3/stack2/classification/schmidhuber/median_filtered/membrane/00_schmidhuber_membrane.tiff';
A = double(imread(filename));
A = A./(max(max(A)));

% C = rgb2gray(A);

% write to another file
% writefile = '/home/thanuja/Dropbox/data/RF_training_edge/I01_trainingLabels.tif';
% writeFilePath = '/home/thanuja/Dropbox/data2/raw';
% writeFilePath = '/home/thanuja/Dropbox/data2/probabilities/neuron';
% writeFilePath = '/home/thanuja/Dropbox/data2/probabilities/neuron';
writeFilePath = '/home/thanuja/Dropbox/data2/probabilities/membrane';
writeFileName = '101.png';

dimx = 250;
dimy = 250;

startRow = 200;
stopRow = startRow -1 + dimy;

startCol = 200;
stopCol = startCol - 1 + dimx;

numDim = 3;

writeFileName = fullfile(writeFilePath,writeFileName);
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
        
