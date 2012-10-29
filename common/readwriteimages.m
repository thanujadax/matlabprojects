% Read png file
% filename = '/home/thanuja/matlabprojects/data/mitoData/labels/stem1_membranes.png';
filename = '/home/thanuja/matlabprojects/data/mitoData/stem1.png';
fmt = 'png';
A = imread(filename, fmt);
C = rgb2gray(A);
% write to another file
writefile = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
B = C(200:248,100:148,1);
imwrite(B,writefile,'png')


