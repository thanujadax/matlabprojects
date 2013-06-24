% load and visualize a set of 3D sections of images
numImg = 20;
dimx = 1024;
dimy = 1024;
imageStack = zeros(dimy,dimx,numImg);
imgDir = '/home/thanuja/Dropbox/data/em_2013january/raw/';
for a = 0:(numImg-1)
   filename = [imgDir num2str(a,'%02d') '.tif'];
   img = imread(filename);
   % add to the 3D image stack
   i = a+1;
   imageStack(:,:,i) = img;
end

yplaneInd = 3;
yplane = zeros(numImg,dimx);
for i=1:numImg
    yplane(i,:) = imageStack(yplaneInd,:,i);
end
figure;imagesc(yplane);colormap('gray')