% input = '/home/thanuja/Dropbox/data/mitoData/stem1_256by256.png';
input = 'stem_256x_t02_V.png';
img = imread(input);
sigma = 0.8;
maskSize = 3;
imSmooth = gaussianFilter(img,sigma,maskSize);
L = watershed(img,8);
figure;imagesc(L);