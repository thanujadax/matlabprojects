% Hough transformation based analysis

%% parameters
imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
houghSupNHood = [5 5];
rhoResolution = 0.5;
thetaRange = -90:0.5:89.5;
houghThresh = 0.5;          % fraction of max(H)
% read image


img = double(imread(imagePath))/255;
figure(1);
imagesc(img);
colormap('gray');
title('original')
% invert
img = invertImage(img);
% rescale 0 - 1
img = img/max(max(img));

figure(2);
imagesc(img);
title('inverted input')
colormap('gray');
% extract edges
BW = edge(img,'canny');

figure(3); 
imagesc(BW);
title('edge detection - canny')
colormap('gray');

% hough transform
[H,T,R] = hough(BW,'RhoResolution',rhoResolution,'Theta',thetaRange);
figure(4)
imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
title('Hough transform');
xtitle();
ytitle();
% get hough peaks
P  = houghpeaks(H,10,'threshold',ceil(0.4*max(H(:))),'NHoodSize',houghSupNHood);
% P = [(row,col), ...]
axis on, axis normal, hold on;
plot(T(P(:,2)),R(P(:,1)),'s','color','white');

