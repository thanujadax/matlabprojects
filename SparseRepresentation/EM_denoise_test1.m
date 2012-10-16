% Thanuja 06.07.2012
% Denoising EM image using K-SVD

clear
%% Parameters
removedEdgeSize = 0; % for preprocessing (removes the dark edges of the section)
bb = 16; % block size
RR = 4; % redundancy factor
% K = RR*bb^2; % number of atoms in the dictionary
%K = 64;
K=400;
maxNumBlocksToTrainOn = 1000;
maxBlocksToConsider = 5000;
sigma = 40;
pathForImages = '/home/thanuja/matlabprojects/data/mitoData/';
% pathForImages = '/home/thanuja/matlabprojects/data/isbi_data/traininglabelspng/';
imageName = 'stem1_256by256.png';
% imageName = 'stemlbl1.png';

%% Get input image
[IMin_0,pp]=imread(strcat([pathForImages,imageName]));

%% preprocessing (to remove the dark edges at the bottom and to the right)
if (1)
    imSize = size(IMin_0);
    resizeDim = [imSize(1)-removedEdgeSize imSize(2)-removedEdgeSize];
    IMin = IMin_0(1:resizeDim(1), 1:resizeDim(2), 1);
end

%% fixing input image format
IMin=im2double(IMin);
if (length(size(IMin))>2)
    IMin = rgb2gray(IMin);
end
if (max(IMin(:))<2)
    IMin = IMin*255;
end

%% Denoising using K-SVD
% dictionary is trained on the noisy image
[IoutAdaptive,output] = denoiseImageKSVD(IMin, sigma,K,bb,'maxNumBlocksToTrainOn',maxNumBlocksToTrainOn,'maxBlocksToConsider',maxBlocksToConsider);

% load('Dictionary3_256')
% [IoutAdaptive,output] = denoiseImageKsvdReuseDict(IMin, sigma,K,bb,Dictionary,'maxNumBlocksToTrainOn',maxNumBlocksToTrainOn);


%PSNROut = 20*log10(255/sqrt(mean((IoutAdaptive(:)-IMin0(:)).^2)));
figure();
imshow(IMin,[]); 
title('Noisy image');
figure();
imshow(IoutAdaptive,[]); 
title('Clean Image by Adaptive dictionary');
figure;
I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
title('The dictionary trained on patches from the noisy image');

