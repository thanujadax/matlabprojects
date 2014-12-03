function image2_patch = getCorrespondingPatch(image2,image1_patch)

% find the corresponding patch image2_patch in image2 which is mostly
% correlated with image1_patch

% Inputs:
%   image2: image in which we should find the patch
%   image1_patch: 

% Output:
%   image2_patch: 

c = normxcorr2(image1_patch,image2);

% offset found by correlation
[max_c, imax] = max(abs(c(:)));
[ypeak, xpeak] = ind2sub(size(c),imax(1));
corr_offset = [(xpeak-size(image1_patch,2))
               (ypeak-size(image1_patch,1))];


% total offset
xoffset = corr_offset(1);
yoffset = corr_offset(2);