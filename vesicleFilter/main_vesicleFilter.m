%% read tiff stack into 3D matrix
% tif file
inputImageVolumeFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/Data/raw/GAS1_mss1_small_50.tif';
outputImageVolumeFileName = '/home/thanuja/Dropbox/PROJECTS/vesicleFilter/small50_1.tif';

imageInfo = iminfo(inputImageVolumeFileName);
%% convolution with 3D structure element (vesicle)
vesicleElement3D = getVesicleElement(innerRadius,outerRadius);

% convn (or fftn for fft based convolution)

% interpn: to interpolate the convolution result at the points you want

%% visualize
imwrite(convolutionOutput,outputImageVolumeFileName,'tif');
% currently to visualize the tiff volume, we are using fiji