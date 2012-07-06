% Thanuja 06.07.2012
% Denoising EM image using K-SVD

clear
%% Parameters
bb = 8; % block size
RR = 4; % redundancy factor
K = RR*bb^2; % number of atoms in the dictionary

sigma = 25;
pathForImages = '/home/thanuja/matlabprojects/data/mitoData';
imageName = '';