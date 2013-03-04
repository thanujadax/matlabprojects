% watershed segmentation
inputfile = '/home/thanuja/Dropbox/data/testImg/stem_256x_t02_V.png';
% inputfile = '/home/thanuja/Dropbox/data/testImg/testMem1_V.png';
I=imread(inputfile);

% I=rgb2gray(I);
thresh = 0.3;
% threshImg = output3(:,:,3);
threshImg = I;
maxResp = max(max(threshImg));
lowInd = threshImg<thresh*maxResp;
highInd = threshImg>=thresh*maxResp;
threshImg(lowInd) = 0;
threshImg(highInd) = 1;
%I = rgb2gray(RGBimg3);
I = threshImg;

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);


%# Normalize.
g = gradmag - min(gradmag(:));
g = g / max(g(:));

th = graythresh(g); %# Otsu's method.
% a = imhmax(g,th/2); %# Conservatively remove local maxima.
a = imhmax(g,th/2);
% th = graythresh(a);
% b = a > th/4; %# Conservative global threshold.
% c = imclose(b,ones(6)); %# Try to close contours.
c = imclose(a,ones(3));
d = imfill(c,'holes'); %# Not a bad segmentation by itself.
%# Use the rough segmentation to define markers.
g2 = imimposemin(g, ~ imdilate( bwperim(a), ones(3) ));
% g2 = imimposemin(g, ~a);
% L = watershed(g2);
L = watershed(a);


figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
% L = watershed(gradmag);
% Lrgb = label2rgb(L);
figure, imagesc(L), title('Watershed transform of gradient magnitude')
figure, imshow(L), title('Watershed transform of gradient magnitude')

% L2(lowInd) = 0;
% figure, imshow(L2), title('Watershed transform of gradient magnitude (Lrgb)')