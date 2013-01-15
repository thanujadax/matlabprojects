% Edge detection using twin-shaded bar templates
%% parameters
displayIntermediateFigures=1;
%imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
% imagePath = 'testImgLines.png';
 
invertImg = 1;      % 1 for membrane images that have to be inverted for Hough transform calculation

grayThresholding = 0;       % 1 if the inverted image should be thresholded
grayThreshold = 0.5;

gaussianFiltering = 0;      % 1 if gaussian filtering should be performed on the input image
sigma = 1;
maskSize = 5;

slidingDist = 1;           % the number of pixels to jump


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% very important tuning parameter
thresholdFraction = 0.5;    % fraction of max(H) to be used as a threshold for peaks
                    % consider the fact that max(H) refers to a global
                    % maximum of H which might overlook smaller line
                    % segments in some patches with less support. 0.5 is
                    % recommended
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

barLength = 9; % should be odd
barWidth = 5; % should be odd

orientations = 0:10:170;    
withBackground = 0;     % plot the detected bars with the original image in the background

%% input preprocessing
imgIn = double(imread(imagePath))/255;
imgIn = imgIn(1:128,129:256);

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
% rescale: -1 to +1
imgInv = 2 * imgInv/max(max(imgInv)) -1;
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
%% Edge detection using twin-shaded bar templates
barVotes3D = getShadedBarVotes(imgInv,barLength,barWidth,orientations,slidingDist);

% Reconstruction of image with bar position, orientation and vote
% information
output = reconstructTwinBars(barVotes3D,orientations,barLength,barWidth);

