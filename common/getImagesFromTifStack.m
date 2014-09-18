% Read a tiff stack of images and write individual images of chosen size

inputFileName = '/home/thanuja/projects/data/FIB_ziqiang/GroupA_Stack1.tif';
outputDirectory = '/home/thanuja/projects/data/FIB_ziqiang/set1/raw';

writeType = 'tif';

zStart = 1;
zStop = 20;

xStart = 1;
xStop = 512;

yStart = 1;
yStop = 512;

% fname = 'my_file_with_lots_of_images.tif';

info = imfinfo(inputFileName);
num_images = numel(info);

k = 0;
for i = zStart:zStop
    A = uint16(imread(inputFileName, i));
    B = A(yStart:yStop,xStart:xStop);
    writeName = sprintf('%02d.%s',k,writeType);
    writeFileName = fullfile(outputDirectory,writeName);
    k = k+1;
    imwrite(B,writeFileName,writeType)
    figure;imagesc(B);
end