% evaluate each cell individually depending on the requirement.

savePath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/raw';
fileName = '03.png';

saveFileName = fullfile(savePath,fileName);
fmt = 'png';
dimx = 256;
dimy = 256;

A = zeros(dimy,dimx);

%% assign colors
startR = 1;
stopR = 128;

startC = 129;
stopC = 256;

color = 0.6;

A(startR:stopR,startC:stopC) = color;


RGBmat = zeros(dimy,dimx,3);
RGBmat(:,:,1) = A;
RGBmat(:,:,2) = A;
RGBmat(:,:,3) = A;

figure;imshow(RGBmat)

%% assign bw
startR = 1;
stopR = 128;

startC = 1;
stopC = 256;

color = 0.9;

A(startR:stopR,startC:stopC) = color;

RGBmat = A;
figure;imshow(RGBmat)
%% save
imwrite(RGBmat,saveFileName,fmt);