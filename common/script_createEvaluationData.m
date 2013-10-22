% inputfilename = '/home/thanuja/Dropbox/data/em_2013january/neurons/06.tiff';
% writeFilePath = '/home/thanuja/Dropbox/data/evaldata2/labels/';

inputfilename = '/home/thanuja/Dropbox/data/em_2013january/raw/06.tif';
writeFilePath = '/home/thanuja/Dropbox/data/evaldata2/input/';

fmt = 'tif';
A = imread(inputfilename, fmt);

dimx = 512;
dimy = 512;

[sizeR,sizeC] = size(A);

numBlocks = sizeR/dimx;

startRow = 1;
stopRow = dimy;

startCol = 1;
stopCol = dimx;

k = 00; % file index
for i=1:numBlocks
    for j=1:numBlocks
        B = A(startRow:stopRow,startCol:stopCol,:);
        % convert RGB into integer labels
        % intLabelsB = im2uint16(B);
        
        % writeName = sprintf('I%02d_neuronLabels06.tif',k);
        
        writeName = sprintf('I%02d_raw06.tif',k);
        
        writeFileName = strcat(writeFilePath,writeName);
        disp(writeFileName)
        % imwrite(intLabelsB,writeFileName,'tif')
        imwrite(B,writeFileName,'tif')
        % figure;imshow(B);
        
        startCol = stopCol + 1;
        stopCol = stopCol + dimx;
        k = k + 1;
    end
    startRow = stopRow + 1;
    stopRow = stopRow + dimy;
    startCol = 1;
    stopCol = dimx;
end