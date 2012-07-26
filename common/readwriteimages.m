% Read png file
filename = '/home/thanuja/matlabprojects/data/mitoData/labels/stem1_membranes.png';
fmt = 'png';
A = imread(filename, fmt);
C = rgb2gray(A);
% write to another file
writefile = '/home/thanuja/matlabprojects/data/mitoData/labels/stem1_memlb_256_2';
B = C(1:256,257:512,1);
imwrite(B,writefile,'png')
