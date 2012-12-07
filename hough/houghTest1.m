% Hough transformation based analysis

%% parameters
imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
houghSupNHood = [5 5];
rhoResolution = 0.5;
thetaRange = -90:0.5:89.5;
houghThresh = 0.5;          % fraction of max(H)
maxLines = 100;
fillGap = 3;                % line segments with a gap less than this will be merged
minLength = 2;              % lines with length below this will be discarded
% read image

%% input 
img = double(imread(imagePath))/255;
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
%% edge detection - canny
% extract edges
BW = edge(imgInv,'canny');

figure(3); 
imagesc(BW);
title('edge detection - canny')
colormap('gray');

%% Hough transform
[H,T,R] = hough(BW,'RhoResolution',rhoResolution,'Theta',thetaRange);
figure(4)
imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
title('Hough transform');
xlabel('theta');
ylabel('rho');
% get hough peaks
P  = houghpeaks(H,maxLines,'threshold',ceil(0.4*max(H(:))),'NHoodSize',houghSupNHood);
% P = [(row,col), ...]

% plot the hough space
axis on, axis normal, hold on;
plot(T(P(:,2)),R(P(:,1)),'s','color','white');
% plot the maxima
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','green');

%% inverse Hough - finding the lines 
lines = houghlines(BW,T,R,P,'FillGap',fillGap,'MinLength',minLength);
figure, imshow(BW), hold on
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



