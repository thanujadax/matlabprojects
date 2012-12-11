% Hough transform based line segment detection with fixed length line
% segments

%% parameters
imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
% imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
rhoResolution = 1;
thetaRange = -90:0.5:89.5;

grayThresholding = 1;       % 1 if the inverted image should be thresholded
grayThreshold = 0.5;

gaussianFiltering = 0;      % 1 if gaussian filtering should be performed on the input image
sigma = 1;
maskSize = 5;

bb = 16;                    % patch size
slidingDist = 8;            % sliding distance for the overlap of patches

%% input preprocessing
imgIn = double(imread(imagePath))/255;
if(size(size(imgIn),2)>2)
    img = imgIn(:,:,1);
else
    img = imgIn;
end
figure(1);
imagesc(img);
colormap('gray');
title('original')
% invert
imgInv = invertImage(img);
% rescale 0 - 1
imgInv = imgInv/max(max(imgInv));

figure(2);
imagesc(imgInv);
title('inverted input')
colormap('gray');

% thresholding
if(grayThresholding == 1)
 imgInv = simpleThreshold(imgInv,grayThreshold);
 figure(8);
 imagesc(imgInv);
 title('inverted input after thresholding')
 colormap('gray');
end

% gaussian smoothening
if(gaussianFiltering==1)
    imgInv = gaussianFilter(imgInv,sigma,maskSize);
    figure(5);
    imagesc(imgInv);
    title('gaussian smoothening');
    colormap('gray');
end

%% Localized Hough transform

% for each image block (overlapping)
% calculate the Hough space and store it in a structure
[localizedHoughSpaces,R,T] = getLocalHoughs(imgInv,rhoResolution,thetaRange,bb,slidingDist);

%% Inverse Hough transform
% Reconstruct the global line sketch using the local Hough spaces
globalHoughLines = localToGlobalHoughLines(localizedHoughSpaces,R,T);
