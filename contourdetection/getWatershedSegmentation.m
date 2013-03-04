function L = getWatershedSegmentation(inputImg,thresh)
% Inputs:
%   inputfilename - path of the input file to be segmented
%   thresh - pixel values above this would be considered as 1 when
%   converted into binary


threshImg = inputImg;
maxResp = max(max(threshImg));

% converting image to binary
lowInd = threshImg<thresh*maxResp;
highInd = threshImg>=thresh*maxResp;
threshImg(lowInd) = 0;
threshImg(highInd) = 1;
I = threshImg;

% calculate gradient
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);


%# Normalize.
g = gradmag - min(gradmag(:));
g = g / max(g(:));

th = graythresh(g); %# Otsu's method for gradient image
% a = imhmax(g,th/2); %# Conservatively remove local maxima.
a = imhmax(g,th/2);
% th = graythresh(a);
L = watershed(a);
