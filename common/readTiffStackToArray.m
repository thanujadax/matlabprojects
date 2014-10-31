function inputImageStack = readTiffStackToArray(inputImageVolumeFileName)

info = imfinfo(inputImageVolumeFileName);
numImages = numel(info);

imHeight = info.Height;
imWidth = info.Width;
inputImageStack = zeros(imHeight,imWidth,numImages);

for k = 1:numImages
    inputImageStack(:,:,k) = imread(inputImageVolumeFileName, k, 'Info', info);
end