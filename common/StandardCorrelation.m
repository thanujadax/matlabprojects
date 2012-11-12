function ima = StandardCorrelation(ima1,ima2)

% Reads in both images
s1 = size(ima1);
s2 = size(ima2);

% Find the max size of the images
s = max(s1,s2);

% Allocates it to the output image
sx = s(1)
sy = s(2)

% Copy original images into bigger images
image1 = zeros(sx,sy);
image1(sx/2-s1(1)/2+1:sx/2+s1(1)/2,sy/2+1-s1(2)/2:sy/2+s1(2)/2) = ima1;
image2 = zeros(sx,sy);
image2(sx/2-s2(1)/2+1:sx/2+s2(1)/2,sy/2+1-s2(2)/2:sy/2+s2(2)/2) = ima2;

% Calculates FFT of bith images
f1 = fft2(fftshift(image1));

f2 = fft2(fftshift(image2));

% Calculates correlation
corr = fftshift(ifft2(f1.*conj(f2)));

% Normalises it with autocorrelation
corr = corr ./ max(max(ifft2(f2.*conj(f2))));
ima = abs(corr);
imagesc(ima);