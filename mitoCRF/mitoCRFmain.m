%% mitoCRFMain.m
% main script
%% Parameters and file paths
FileSTEM = '../data/mitoData/stem1.tiff';
FileTEM = '../data/mitoData/abd.tiff';
FileISBI = '../data/isbi_data/train-labels.tif';
OutputTiff = './outputs/probmap_schmidhuber.tiff';


%% Read Tiff image stack
TifImgStack = readTiffStack(FileISBI);

%% Invert 
% newmat2 = uint8(ones(size(TifImgStack,1),size(TifImgStack,2)) .* 255);
% invertedImage = newmat2 - TifImgStack(:,:,1) ;

%% Graph Cut
% outputGC = CMF_Cut(TifImgStack(:,:,1),0.2);
% imshow(outputGC);

%% Median filtering
% medfil = medfilt2(TifImgStack(:,:,1), [3 3]);
% imshow(medfil);

%% Write Tiff image stack
% % writeTiffStack(TifImgStack, OutputTiff);
% writeTiffStack(invertedImage, OutputTiff);

%% Write png image
A = TifImgStack(1:512,1:512,1);
pngfile ='../data/isbi_data/traininglabelspng/stemlbl1.png';
imwrite(A,pngfile,'png')

%% show img
% inImage = TifImgStack(:,:,1);
% imSize = size(inImage);
% 
% removedEdgeSize = 4;
% resizeDim = [imSize(1)-removedEdgeSize imSize(2)-removedEdgeSize];
% inImage2 = inImage(1:resizeDim(1), 1:resizeDim(2));
% figure();
% imshow(inImage);
% figure();
% imshow(inImage2);