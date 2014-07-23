function outputImage = fillHoles(inputImage)

outputImage = inputImage;

if (size(inputImage,3)==3)
    im_1 = inputImage(:,:,1);
    im_2 = inputImage(:,:,2);
    im_3 = inputImage(:,:,3);
    
    outputImage(:,:,1) = imfill(im_1);
    outputImage(:,:,2) = imfill(im_2);
    outputImage(:,:,3) = imfill(im_3);
    
elseif(size(inputImage,3)==1)
    outputImage = imfill(inputImage);
    
else
    disp('Error fillHoles.m. Check dimensions of input image')
end