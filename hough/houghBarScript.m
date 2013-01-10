% Hough: detecting oriented bars

%% parameters
displayIntermediateFigures=0;
% imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
% imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
imagePath = 'testImgLines.png';
 
invertImg = 0;      % 1 for membrane images that have to be inverted for Hough transform calculation

rhoResolution = 0.5;
thetaRange = -90:0.5:89.5;
%thetaRange = 0:0.5:90.0;

grayThresholding = 1;       % 1 if the inverted image should be thresholded
grayThreshold = 0.5;

gaussianFiltering = 0;      % 1 if gaussian filtering should be performed on the input image
sigma = 1;
maskSize = 5;

bb = 32;                    % patch size
slidingDist = 1;           % the number of pixels to jump

maxLinesPerPatch = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% very important tuning parameter
thresholdFraction = 0.8;    % fraction of max(H) to be used as a threshold for peaks
                    % consider the fact that max(H) refers to a global
                    % maximum of H which might overlook smaller line
                    % segments in some patches with less support. 0.5 is
                    % recommended
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
houghSupNHood = [5 5];      % suppression neighborhood at each identified peak
fillGap = 2;                % fill gaps smaller than this to combine two collinear lines    
minLength = 20;              % minimum length of lines to be detected

smoothenH = 1;      % if each local H should be smoothened using a gaussian filter
sigmaH = 0.3;
maskSizeH = 3;

barLength = 9;
barWidth = 3;
orientations = [0 45 90 135];    % can either be 4 or 8

withBackground = 0;     % plot the detected bars with the original image in the background

%% input preprocessing
imgIn = double(imread(imagePath))/255;
% imgIn = imgIn(1:128,1:128);

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
%% Hough
% perform Hough type processing (voting) for oriented bars
houghSpace3D = houghBars(imgInv,barLength,barWidth,orientations,slidingDist);
% houghSpace3D [row col orientation]

% peak detection
peaks3D = houghBarPeaks(houghSpace3D,orientations,thresholdFraction...
                            ,slidingDist,barLength,barWidth);  

% draw the detected bars on the image
output = reconstructHoughBars(peaks3D,orientations,barLength,barWidth);
figure(101);imagesc(output);title('reconstruction using oriented bars')

