% Thanuja 26.06.2012
%% input data
infile1 = '/home/thanuja/matlabprojects/data/curves2/curve2.jpg';
infile2 = '/home/thanuja/matlabprojects/data/mitoData/stem1.tiff';
infile3 = '/home/thanuja/matlabprojects/data/curves/curves_noise.png';
%img1 = imread(infile1,'jpg');
% img1 = rgb2gray(imread(infile3, 'png'));
img1 = imread(infile2, 'tif');

img1 = double(img1)./255;

%% filter bank parameters
max_freq = 0.01;
num_freq = 1;
num_orientation = 4;

%% create filter bank
bank=sg_createfilterbank(size(img1), max_freq , num_freq, num_orientation,'verbose',1);
%% filtering and output generation
% filter image with filter bank
r=sg_filterwithbank(img1,bank,'method',1);

% convert resulting structure r into 3x3 matrix
% dimensions of m = (resolution_X) x (resolution_Y) x (# gabor filters)
m = sg_resp2samplematrix(r);

% sum up all the responses and view
FIG1 = figure();
imagesc(abs(sum(m,3))); colormap(FIG1,gray);







