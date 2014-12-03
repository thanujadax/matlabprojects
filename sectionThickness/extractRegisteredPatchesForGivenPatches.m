function image2_patchesRegistered = extractRegisteredPatchesForGivenPatches(...
                                    image2,image1_patches)
                                
% Inputs:
%   image2 - image to be used to make patches
%   image1_patches - image patches to be registered with

% Outputs:

[patchSizeR,patchSizeC,numPatches] = size(image1_patches);
image2_patchesRegistered = zeros(patchSizeR,patchSizeC,numPatches);

for i=1:numPatches
    image1_patch_i = image1_patches(:,:,i);
    image2_patch_i = getCorrespondingPatch(image2,image1_patch_i);
    image2_patchesRegistered(:,:,i) = image2_patch_i;
end


