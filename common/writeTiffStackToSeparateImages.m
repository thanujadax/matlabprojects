function writeTiffStackToSeparateImages(inputImageVolumeFileName,saveFilePath,fileType)

inputImageStack = readTiffStackToArray(inputImageVolumeFileName);

numImg = size(inputImageStack,3);

%fileType = 'tif';

for i=1:numImg
    imageName = sprintf('%d.%s',(i-1),fileType);
    outputFullFile = fullfile(saveFilePath,imageName);
    disp(outputFullFile);
    image_i = inputImageStack(:,:,i);
    imagesc(image_i);
    imwrite(image_i,outputFullFile);
end