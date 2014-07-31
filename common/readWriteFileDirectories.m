% Read files from directory and write processed files to another directory

%% file paths and names

% filename = '/home/thanuja/projects/drosophila-l3/stack2/raw/00.tif';
% filename = '/home/thanuja/projects/drosophila-l3/stack2/classification/schmidhuber/median_filtered/neurons/neurons0000.png';
% filename = '/home/thanuja/projects/drosophila-l3/stack2/classification/schmidhuber/median_filtered/membrane/00_schmidhuber_membrane.tiff';
% filename = '/home/thanuja/projects/drosophila-l3/stack2/groundtruth/result_0000.tiff';

% writeFilePath = '/home/thanuja/projects/toyData/set9/groundtruth';
% writeFilePath = '/home/thanuja/projects/toyData/set9/membranes';
% writeFilePath = '/home/thanuja/projects/toyData/set9/neurons';
% writeFilePath = '/home/thanuja/projects/toyData/set9/raw';

% inputDirectoryName = '/home/thanuja/projects/drosophila-l3/stack2/raw/';
% inputDirectoryName =  '/home/thanuja/projects/drosophila-l3/stack2/classification/schmidhuber/median_filtered/neurons/';
% inputDirectoryName = '/home/thanuja/projects/drosophila-l3/stack2/classification/schmidhuber/median_filtered/membrane/';
inputDirectoryName = '/home/thanuja/projects/drosophila-l3/stack2/groundtruth/';

imgFileExtension = 'tiff';

% writeFilePath = '/home/thanuja/projects/inputData/trainingHalf/raw/';
% writeFilePath = '/home/thanuja/projects/inputData/trainingHalf/neurons/';
% writeFilePath = '/home/thanuja/projects/inputData/trainingHalf/membranes/';

% writeFilePath = '/home/thanuja/projects/inputData/trainingHalf';
writeFilePath = '/home/thanuja/projects/inputData/testingHalf';

% writeSubDir = 'raw';
% writeSubDir = 'neurons';
% writeSubDir = 'membranes';
writeSubDir = 'groundtruth';

writeType = 'png';

%% param
dimx = 1024;
dimy = 512;

startRow = 513;
stopRow = startRow -1 + dimy;

startCol = 1;
stopCol = startCol - 1 + dimx;

numDim = 3;
%% read images from input directory
imageFileDirectoryPath = strcat(inputDirectoryName,'*');
imageFileDirectoryPath = strcat(imageFileDirectoryPath,imgFileExtension);
imageFiles = dir(imageFileDirectoryPath);
nfiles = length(imageFiles);    % Number of files found

%% write to output directory
% A = double(imread(filename));
% A = A./(max(max(A)));
% A = A./255;



k = 00; % file index
for i=1:nfiles
    currentInputFileName = imageFiles(i).name;
    inputFullFile = fullfile(inputDirectoryName,currentInputFileName);
    disp('Input file: ');
    disp(inputFullFile);
    A = double(imread(inputFullFile));
    A = A./255;
    B = A(startRow:stopRow,startCol:stopCol,:);
    writeName = sprintf('%02d.%s',k,writeType);
    writeFileName = fullfile(writeFilePath,writeSubDir,writeName);
    disp('Output file: ')
    disp(writeFileName)
    imwrite(B,writeFileName,writeType)
    figure;imshow(B);
    k = k + 1;
end
        