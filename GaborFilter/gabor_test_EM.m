clear;
clf;

% % Load the image
% I=rgb2gray(imread('lena.jpg','jpg'));

%inputfile = '/home/thanuja/matlabprojects/data/mitoData/schmidhuber1.tif';
inputfile = '/home/thanuja/matlabprojects/data/mitoData/stem1.tiff';

I0=imread(inputfile,'tiff');
I0 = uint8((double(I0)-255) * -1);

load colormaps.mat

% Show the grayscale image
figure(10);
colormap(grayscale);
imshow(I0);
title('original image')

%% parameters
variance = 0.1;
freqAmplitude = 0.001;
orientation = 0;  % pi/2
phase = 0;
linearthresh = 50; % thresholding the Gabor output (magnitude) for denoising
%% Filter the image
[G,GABOUT]=gaborfilter(I0,variance,freqAmplitude,orientation,phase);
%[G,GABOUT]=gaborfilter(I,0.05,0.025,0,0);

% clear I0;

R=real(GABOUT);
I=imag(GABOUT);
M=abs(GABOUT);
P=angle(GABOUT);

clear GABOUT;
%% plots
% % show gabor filter response
% figure(14);
% colormap(jet);
% surf(real(G));
% title('real part of Gabor filter');
% 
% figure(15);
% colormap(jet);
% surf(abs(G));
% title('abs of Gabor filter');

% % Show the filter's outputs
% figure(11);
% colormap(redgreen);
% subplot(2,2,1);
% k=127.5/max(max(abs(R)));
% image(uint8(k*R+127.5));
% subplot(2,2,2);
% k=127.5/max(max(abs(I)));
% image(uint8(k*I+127.5));

% % Show the kernels
% colormap(redgreen);
% subplot(2,2,3);
% image(uint8(127.5*real(G)+127.5));
% subplot(2,2,4);
% image(uint8(127.5*imag(G)+127.5));

% Show the magnitudes
figure();
colormap(grayscale);
k=255/max(max(M));
absgaborout = k*M;
image(uint8(absgaborout));
title('abs of Gabor output (magnitude)');

% % Show the phases
% figure(13);
% colormap(redgreen);
% k=127.5/pi;
% image(uint8(k*P+127.5));
% title('phase of Gabor output');

%% post processing
% thresholding Gabor output (magnitude)
gabormag = absgaborout;
smallind = find(absgaborout<70);
gabormag(smallind)=0;
figure();
colormap(gray);
%imshow(gabormag);
image(gabormag);
title('thresholding gabor magnitude 70')

% thresholding Gabor output (magnitude)
%gabormag = absgaborout;
whitearea = find(absgaborout>70);
newimage = zeros(size(absgaborout));
newimage(whitearea)=255;
figure();
colormap(gray);
%imshow(gabormag);
image(newimage);
title('thresholding gabor magnitude 70 (whitearea)')
