% Hough transformation based analysis

%% parameters
displayIntermediateFigures = 0;
imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
% imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
houghSupNHood = [5 5];
rhoResolution = 1;
thetaRange = -90:0.5:89.5;

houghThresh = 0.5;          % fraction of max(H)
maxLines = 100;
fillGap = 2;                % line segments with a gap less than this will be merged
minLength = 4;              % lines with length below this will be discarded

gaussianFiltering = 0;      % 1 if gaussian filtering should be performed on the input image
sigma = 1;
maskSize = 5;

grayThresholding = 1;       % 1 if the inverted image should be thresholded
grayThreshold = 0.5;

edgeDetection = 0;          % 1 if edge detection should be performed

gaussSmoothH = 0;           % 1 if H should be smoothed by Gaussian
sigmaH = 2;
maskSizeH = 5;
% read image

%% input 
imgIn1 = double(imread(imagePath))/255;
imgIn = imgIn1(1:24,1:16);
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
imgInv = invertImage(img);
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
%% gaussian smoothening
if(gaussianFiltering==1)
    if(displayIntermediateFigures)
        imgInv = gaussianFilter(imgInv,sigma,maskSize);
        figure(5);
        imagesc(imgInv);
        title('gaussian smoothening');
        colormap('gray');
    end
end

%% edge detection - canny
% extract edges
if(edgeDetection==0)
    BW = imgInv;
else
    BW = edge(imgInv,'canny');
    if(displayIntermediateFigures)
        figure(3); 
        imagesc(BW);
        title('edge detection - canny')
        colormap('gray');
    end
end

%% Hough transform
% [H,T,R] = hough(BW,'RhoResolution',rhoResolution,'Theta',thetaRange);
[H,T,R] = hough2(BW,rhoResolution,thetaRange);
if(displayIntermediateFigures)
    figure(4)
    imshow(H,[],'XData',T,'YData',R,...
                'InitialMagnification','fit');
    title('Hough transform');
    xlabel('theta');
    ylabel('rho');
end

%% gaussian smoothening H
if(gaussSmoothH==1)
    H = gaussianFilter(H,sigmaH,maskSizeH);
    H = ceil(H);
    imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
    title('Hough transform after smoothening');
    xlabel('theta');
    ylabel('rho');
end

% get hough peaks
P  = houghpeaks(H,maxLines,'threshold',ceil(0.4*max(max(H)),'NHoodSize',houghSupNHood);
% P = [(row,col), ...]
if(displayIntermediateFigures)
    % plot the hough space
    axis on, axis normal, hold on;
    plot(T(P(:,2)),R(P(:,1)),'s','color','white');
    % plot the maxima
    x = T(P(:,2)); y = R(P(:,1));
    plot(x,y,'s','color','green');
end

%% inverse Hough - finding the lines 
lines = houghlines(BW,T,R,P,'FillGap',fillGap,'MinLength',minLength);
figure(7), imagesc(BW), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

% highlight the longest line segment
% plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');



