% Hough: detecting oriented bars

%% parameters
displayIntermediateFigures=1;
%imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
% imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
%imagePath = '/home/thanuja/Dropbox/data/em_2013january/samples/raw00_256.png';
imagePath = '/home/thanuja/Dropbox/data/mitoData/stem1_512.png';
% imagePath = '/home/thanuja/Dropbox/data/em_2013january/samples/raw00_512.png';
%imagePath = 'testImgLines.png';
% imagePath = 'testImgLines3.png';
% imagePath = 'testImgGauss.png';
%imagePath = 'testCirc.png';
%imagePath = 'testImgGauss.png';
%imagePath = 'testCirc.png';

Hthresh = 0.4; % pixels above this value will be used for hough voting
 
invertImg = 1;      % 1 for membrane images that have to be inverted for Hough transform calculation

grayThresholding = 0;       % 1 if the inverted image should be thresholded
grayThreshold = 0.5;

gaussianFiltering = 0;      % 1 if gaussian filtering should be performed on the input image
sigma = 1;
maskSize = 5;

slidingDist = 1;           % the number of pixels to jump

lineWidth = 1;

threshFrac = 0.3;

% for Gaussian kernel
% barLength = 23; % should be odd
% barWidth = 7; % should be odd

% for bars
barLength = 11; % should be odd
barWidth = 6; % 
negLines = 3; % number of negative lines per side

orientations = 0:10:170;    
withBackground = 0;     % plot the detected bars with the original image in the background

sigmaDeriv = 0.5;   % for the gaussian derivative (to produce edge map)

% for gaussian kernel
sigX = 40;
sigY = 6;

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
    imshow(img);
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
    imshow(imgInv);
    title('inverted input')
    colormap('gray');
end
% thresholding
if(grayThresholding == 1)
 imgInv = simpleThreshold(imgInv,grayThreshold);
 if(displayIntermediateFigures)
     figure(8);
     imshow(imgInv);
     title('inverted input after thresholding')
     colormap('gray');
 end
end

% gaussian smoothening
if(gaussianFiltering==1)
    imgInv = gaussianFilter(imgInv,sigma,maskSize);
    if(displayIntermediateFigures)
        figure(5);
        imshow(imgInv);
        title('gaussian smoothening');
        colormap('gray');
    end
end
%% convolution
% convolution using oriented gaussian kernels
display('Computing 3D orientation score space...');
t0 = cputime;
% orientedScoreSpace3D = convolveOrientedGauss_P(imgInv,barLength,barWidth,...
%             orientations,sigX,sigY);
orientedScoreSpace3D = convolveOrientedBars_P(imgInv,barLength,barWidth,...
            orientations,negLines);
t1 = cputime;
display('3D orientation score space computed!');
dt = t1 - t0;
str = sprintf('Time taken for score space creation = %0.5f s',dt);
disp(str);


%%
[output3 RGBimg3] = reconstructHSVgauss_mv(orientedScoreSpace3D,orientations,barLength,barWidth,threshFrac);
titlestr = sprintf('threshold percentage = %f',threshFrac);
figure;imshow(RGBimg3);title(titlestr)
% batch processing
savefilepath = '/home/thanuja/Dropbox/RESULTS/orientations/thresholding2/';
for i=0:5:99
    % run reconstruction for threshold = i/100
    threshFrac = i/100;
    [output3 RGBimg3] = reconstructHSVgauss_mv(orientedScoreSpace3D,orientations,barLength,barWidth,threshFrac);
    % save it in a folder
    savefilename = sprintf('threshPercent%d.png',i);
    savefilename = strcat(savefilepath,savefilename);
    titlestr = sprintf('threshold = %f',threshFrac);
    h = figure;imshow(RGBimg3);title(titlestr);    
    saveas(h,savefilename);
end


% [output RGBimg] = reconstructHSVbars(orientedScoreSpace3D,orientations,barLength,barWidth,threshFrac);
% writeFile1 = '/home/thanuja/Dropbox/RESULTS/hough/orientations/reconst_hough256raw00_L11_1.png';
% imwrite(RGBimg,writeFile1,'png');
% [output2 RGBimg2] = reconstructHSVlines(houghSpace3D,orientations,barLength,lineWidth,threshFrac);
% writeFile2 = '/home/thanuja/Dropbox/RESULTS/hough/orientations/reconst_hough256raw00_L11_lines_1.png';
% imwrite(RGBimg,writeFile2,'png');
