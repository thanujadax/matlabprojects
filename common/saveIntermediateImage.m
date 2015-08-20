function saveIntermediateImage(ImgIn,rawImageID,intermediateImgDescription,...
    saveIntermediateImagesPath)

outputFileName = sprintf('%s_%s.png',rawImageID,intermediateImgDescription);
outputFullFile = fullfile(saveIntermediateImagesPath,outputFileName);
imwrite(ImgIn,outputFullFile)