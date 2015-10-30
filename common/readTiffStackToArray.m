function inputImageStack = readTiffStackToArray(inputImageVolumeFileName)

info = imfinfo(inputImageVolumeFileName);
numImages = numel(info);

imHeight = info.Height;
imWidth = info.Width;
inputImageStack = zeros(imHeight,imWidth,numImages);

if(strcmp(info(1).Format,'png'))

    tmp1 = imread(inputImageVolumeFileName);
    
    if(size(tmp1,3)==3)
        inputImageStack(:,:) = rgb2gray(tmp1);
    else
        inputImageStack(:,:) = tmp1;
    end
    
else

    for k = 1:numImages
        inputImageStack(:,:,k) = imread(inputImageVolumeFileName, k, 'Info', info);
    end
    
end