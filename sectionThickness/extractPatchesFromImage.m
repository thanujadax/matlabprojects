function imagePatches = extractPatchesFromImage(...
    inputImageName,patchSizeC,patchSizeR,maxNumPatches,overlap)

% extract non overlapping patches

% TODO: extend to have overlapping patches as well

inputImage = double(imread(inputImageName));
[sizeR,sizeC] = size(inputImage);
if(patchSizeC>sizeC || patchSizeR>sizeR)
    disp('Error extractPatchesFromImage: patch size exceeds image size!')
    imagePatches = -1;
else
    rStartVect = 1:patchSizeR:sizeR;
    cStartVect = 1:patchSizeC:sizeC;
    
    numPatches_1 = numel(rStartVect) * numel(cStartVect);
    imagePatches = zeros(patchSizeR,patchSizeC,numPatches_1);
    
    for i=1:numPatches_1
        rStart = rStartVect(i);
        rStop = rStart + patchSizeR -1;
        cStart = cStartVect(i);
        cStop = cStart + patchSizeC -1;
        imagePatches(:,:,i) = inputImage...
            (rStart:rStop,cStart:cStop);
    end
% if numPatches > maxNumPatches, return a random subset

% if numPatches < maxNumPatches, extract overlapping patches


end
