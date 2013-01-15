% canny main script

image = '/home/thanuja/matlabprojects/data/mitoData/stem1.tiff';
sizeof_guassian_kernel = [3 3];
sigma = 0.5;
thresh = 0.1;
canny(image,sizeof_guassian_kernel,sigma,thresh)