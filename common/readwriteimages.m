% Read png file
% filename = '/home/thanuja/matlabprojects/data/mitoData/labels/stem1_membranes.png';
% filename = '/home/thanuja/matlabprojects/data/mitoData/stem.tiff';
filename = '/home/thanuja/Dropbox/data/em_2013january/raw/00.tif';
fmt = 'tif';
A = imread(filename, fmt);
% C = rgb2gray(A);
% write to another file
% writefile = '/home/thanuja/Dropbox/data/em_2013january/samples/raw00_512.png';
writefile = '/home/thanuja/Dropbox/data/mitoData/emJ_00_350x.png';
B = A(101:450,101:450,1);
imwrite(B,writefile,'png')
figure;imshow(B);
