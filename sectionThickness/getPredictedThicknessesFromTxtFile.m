function predictedThickness = getPredictedThicknessesFromTxtFile(inputFilePath)

% Reads the .txt files saved in inputFilePath. Each of those files should
% correspond to section thickness estimates.

% the output is a m by n matrix where m is the number of different
% estimates corresponding to different txt files. n is the number of
% sections.

% inputFilePath = '/home/thanuja/projects/tests/thickness/similarityCurves/s108';
txtFileNames = '*.txt';
txtFileNames = fullfile(inputFilePath,txtFileNames);

txtFileDir = dir(txtFileNames);

numFiles = length(txtFileDir);

if(numFiles==0)
    error('No txt files found in %s', inputFilePath);
end

% load one txt file to get the number of data points
txtFileName = fullfile(inputFilePath,txtFileDir(1).name);
fileID = fopen(txtFileName,'r');
formatSpec = '%f';
thickness = fscanf(fileID,formatSpec);
fclose(fileID);
numDataPoints = size(thickness,1);

predictedThickness = zeros(numFiles,numDataPoints);

for i=1:numFiles
    txtFileName = fullfile(inputFilePath,txtFileDir(i).name);
    fileID = fopen(txtFileName,'r');
    formatSpec = '%f';
    predictedThickness(i,:) = fscanf(fileID,formatSpec);
    fclose(fileID);
end