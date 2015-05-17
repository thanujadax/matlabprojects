%% Parameters
sigmaGaussianBlur = 2;
maskSizeGaussianBlur = 7;
threshold_pixelIntensity = 0.15; % to generate binary images for myelin labels


%% Import image
inputImageFileName = '/home/thanuja/projects/data/myelin/RFCprobabilities_ilastik/ssSEM/s909/myelin_01.png';
inputImage = imread(inputImageFileName);
figure;imagesc(inputImage);title('Input Image')

%% Gaussian blur
   
gbImage = gaussianFilter(inputImage,sigmaGaussianBlur,maskSizeGaussianBlur);
figure;imagesc(inputImage);title('Gaussian blurred image')

%% Thresholding

thresholdedGbImage = im2bw(gbImage,threshold_pixelIntensity);
figure;imshow(thresholdedGbImage);title('thresholded gaussian blurred image')

%% Get connected components
CC = bwconncomp(thresholdedGbImage);
numPixels = cellfun(@numel,CC.PixelIdxList);
% [biggest,idx] = max(numPixels);
[sortedNumPixelList,cellIDsForSortedPixelsList] = sort(numPixels,'descend');

% pick the largest
% to access the pixel listCC.PixelIdxList{idx}
% biggest CC
pixIndLogical_largest_1 = CC.PixelIdxList{cellIDsForSortedPixelsList(1)};
canvas1 = zeros(size(gbImage));
canvas1(pixIndLogical_largest_1) = 1;
figure;imshow(canvas1);title('largest CC')

% 2nd largest
pixIndLogical_largest_2 = CC.PixelIdxList{cellIDsForSortedPixelsList(2)};
canvas2 = zeros(size(gbImage));
canvas2(pixIndLogical_largest_2) = 1;
figure;imshow(canvas2);title('2nd largest CC')

% 3rd largest
pixIndLogical_largest_3 = CC.PixelIdxList{cellIDsForSortedPixelsList(3)};
canvas3 = zeros(size(gbImage));
canvas3(pixIndLogical_largest_3) = 1;
figure;imshow(canvas3);title('3rd largest CC')

%% Erosion and dilation
SE = strel('square',5);

erodedImage = imerode(gbImage,SE);
figure;imagesc(erodedImage);title('eroded')

dilatedImage = imdilate(gbImage,SE);
figure;imagesc(dilatedImage);title('dilated')

erosionDilation = imdilate(erodedImage,SE);
figure;imagesc(erosionDilation);title('eroded dilated')

dilationErosion = imerode(dilatedImage,SE);
figure;imagesc(erosionDilation);title('dilated eroded')

%% connected component analysis
threshold_pixelIntensity = 0.15;
binaryImage = im2bw(dilationErosion,threshold_pixelIntensity);
figure;imshow(binaryImage);title('binary image after dilation and erosion')

%% watershed
% invertedGbImage = invertImage(gbImage);
% ws0 = watershed(invertedGbImage);
% figure;imagesc(ws0);title('original ws')
% 
% ws = assignRandomIndicesToWatershedTransform(ws0);
% figure;imagesc(ws);title('randomized ws indices')
% 
% watershedColoredEdges = overlayWatershedOnImage(inputImage,ws);
%% threshold

