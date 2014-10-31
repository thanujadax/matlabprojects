function writeTiffStackToSeparateImages(inputImageVolumeFileName)

inputImageStack = readTiffStackToArray(inputImageVolumeFileName);

numImg = size(inputImageStack,3);

fileType = 'tif';

for i=1:numImg
    imageName = sprintf('%d.%s',(i-1),fileType);
end