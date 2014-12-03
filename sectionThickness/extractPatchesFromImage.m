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
    
    rStartVect(end) = [];
    cStartVect(end) = [];
    
    numPatches_1 = numel(rStartVect) * numel(cStartVect);
    imagePatches = zeros(patchSizeR,patchSizeC,numPatches_1);
    
    k = 1;
    
    for i=1:numel(rStartVect)
        rStart = rStartVect(i);
        rStop = rStart + patchSizeR -1;
        for j=1:numel(cStartVect)
            cStart = cStartVect(j);
            cStop = cStart + patchSizeC -1;
            imagePatches(:,:,k) = inputImage...
            (rStart:rStop,cStart:cStop);
            k = k+1;
        end
    end
% if numPatches > maxNumPatches, return a random subset

% if numPatches < maxNumPatches, extract overlapping patches


end
