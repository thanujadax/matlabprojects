function inputImageStack = readTiffStackToArray(inputImageVolumeFileName)

info = imfinfo(inputImageVolumeFileName);
numImages = numel(info);

imHeight = info.Height;
imWidth = info.Width;
inputImageStack = zeros(imHeight,imWidth,numImages);

if(strcmp(info(1).Format,'png'))

    inputImageStack(:,:) = imread(inputImageVolumeFileName);
    
else

    for k = 1:numImages
        inputImageStack(:,:,k) = imread(inputImageVolumeFileName, k, 'Info', info);
    end
    
end