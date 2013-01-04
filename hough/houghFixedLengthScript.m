% Hough transform based line segment detection with fixed length line
% segments

%% parameters
displayIntermediateFigures=0;
% imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
% imagePath = 'testImgLines.png';
 
invertImg = 1;      % 1 for membrane images that have to be inverted for Hough transform calculation

rhoResolution = 0.5;
thetaRange = -90:0.5:89.5;
%thetaRange = 0:0.5:90.0;

grayThresholding = 1;       % 1 if the inverted image should be thresholded
grayThreshold = 0.5;

gaussianFiltering = 0;      % 1 if gaussian filtering should be performed on the input image
sigma = 1;
maskSize = 5;

bb = 32;                    % patch size
slidingDist = 16;           % the number of overlapping pixels

maxLinesPerPatch = 20;
thresholdFraction = 0.5;    % fraction of max(H) to be used as a threshold for peaks
                    % consider the fact that max(H) refers to a global
                    % maximum of H which might overlook smaller line
                    % segments in some patches with less support. 0.5 is
                    % recommended
                    
houghSupNHood = [5 5];      % suppression neighborhood at each identified peak
fillGap = 2;                % fill gaps smaller than this to combine two collinear lines    
minLength = 20;              % minimum length of lines to be detected

smoothenH = 1;      % if each local H should be smoothened using a gaussian filter
sigmaH = 0.3;
maskSizeH = 3;

%% input preprocessing
imgIn = double(imread(imagePath))/255;
imgIn = imgIn(1:128,1:128);

if(size(size(imgIn),2)>2)
    img = imgIn(:,:,1);
else
    img = imgIn;
end
if(displayIntermediateFigures)
    figure(1);
    imagesc(img);
    colormap('gray');
    title('original')
end
% invert
if(invertImg)
    imgInv = invertImage(img);
else
    imgInv = img;
end
% rescale 0 - 1
imgInv = imgInv/max(max(imgInv));
if(displayIntermediateFigures)
    figure(2);
    imagesc(imgInv);
    title('inverted input')
    colormap('gray');
end
% thresholding
if(grayThresholding == 1)
 imgInv = simpleThreshold(imgInv,grayThreshold);
 if(displayIntermediateFigures)
     figure(8);
     imagesc(imgInv);
     title('inverted input after thresholding')
     colormap('gray');
 end
end

% gaussian smoothening
if(gaussianFiltering==1)
    imgInv = gaussianFilter(imgInv,sigma,maskSize);
    if(displayIntermediateFigures)
        figure(5);
        imagesc(imgInv);
        title('gaussian smoothening');
        colormap('gray');
    end
end

%% Localized Hough transform

% for each image block (overlapping)
% calculate the Hough space and store it in a structure
[localHoughSpaces,patchLocations,R,T,maxHoughPeak] ...
    = getLocalHoughSpaces(imgInv,rhoResolution,thetaRange,bb,slidingDist...
            ,smoothenH,sigmaH,maskSizeH);

% Plot the lines from localHoughSpaces
plotLocalHoughSpaces(imgInv,localHoughSpaces,T,R,maxLinesPerPatch,...
                    thresholdFraction,houghSupNHood,fillGap,minLength,...
                    bb,slidingDist);
%% Inverse Hough transform
% Reconstruct the global line sketch using the local Hough spaces
% maxHoughPeak is already detected globally from the previous step and is
% used as a reference for peak detection in individual patches
patchLines = localHoughLines(localHoughSpaces,R,T,imgInv,bb,maxLinesPerPatch,...
            thresholdFraction,houghSupNHood,fillGap,minLength,maxHoughPeak,slidingDist);
% patchLines is a cell array
        
%% Reconstruct local patches from the patch lines
%localPatches = reconstructLocalPatches(patchLines,bb);

%% Reconstruct the image
patchLinesGlobal = patchLinesToGlobal(patchLines,slidingDist,bb);
h = plotHoughLines(patchLinesGlobal,imgInv);


