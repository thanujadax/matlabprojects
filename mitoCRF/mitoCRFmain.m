%% mitoCRFMain.m
% main script
%% Parameters and file paths
FileSTEM = '../data/mitoData/stem.tiff';
FileTEM = '../data/mitoData/abd.tiff';
OutputTiff = './outputs/output.tiff';

%% Read Tiff image stack
TifImgStack = readTiffStack(FileSTEM);

%% Write Tiff image stack
writeTiffStack(TifImgStack, OutputTiff);