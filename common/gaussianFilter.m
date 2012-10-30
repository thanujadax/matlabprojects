function imageOut = gaussianFilter(imageIn,sigma,maskSize)
% returns gaussian filtered image


%# Create the gaussian filter with hsize = [5 5] and sigma = 2
G = fspecial('gaussian',[maskSize maskSize],sigma);
%# Filter it
imageOut = imfilter(imageIn,G,'same');
%# Display
imshow(imageOut)