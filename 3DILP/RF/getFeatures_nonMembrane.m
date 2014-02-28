function [fm]  = getFeatures_nonMembrane(imIn, oriFiltLen, halfWidth_strucEl, csHist)

% Inputs:
%   imIn - input image for which features will be computed. Neighboring
%   sections are NOT considered at this point.
%   oriFiltLen - length of the orientation filter
%   halfwidth_strucEl - halfwidth -1 of the orientation filter
%   csHist

% typical values for the input arguments 
% oriFiltLen = 29;
% halfwidthOfStrElement = 3; % actual width = halfwidth + 1
% csHist = maskSize;

% Output:
%   fm - feature matrix (3D)
%       Each entry in the 3rd dim stands for a a feature
%       Does not take into account the sections before and after.

%% Init

numFeatures = 89; % TODO
fm = single(zeros(size(imIn,1), size(imIn,2), numFeatures));

%% Oriented filter responses
%normalize image contrast
imIn = norm01(imIn);

oriFiltMask = zeros(oriFiltLen,oriFiltLen);
oriFiltMask_center = round(oriFiltLen / 2);
oriFiltMask...
(:,oriFiltMask_center-halfWidth_strucEl:oriFiltMask_center+halfWidth_strucEl)...
= 1;
oriFiltMask = single(oriFiltMask);


fm(:,:,1) = imIn; % this section
% fm(:,:,2) = imIn_1; % section before
% fm(:,:,3) = imIn_2; % section after

imIn = adapthisteq(imIn);

% oriented filter response for this section
rot = filterImageWithMembraneTemplateRotated(imIn, oriFiltMask);

imIn = single(imIn);

fm(:,:,2) = rot(:,:,1);
fm(:,:,3) = rot(:,:,2);
fm(:,:,4) = rot(:,:,3);
fm(:,:,5) = rot(:,:,4);
fm(:,:,6) = rot(:,:,5);
fm(:,:,7) = rot(:,:,6);
fm(:,:,8) = rot(:,:,7);
fm(:,:,9) = rot(:,:,8);

rot = shiftdim(rot,2);

rotMin = min(rot);
fm(:,:,10) = shiftdim(rotMin,1);
clear rotMin;

rotMax = max(rot);
fm(:,:,11) = shiftdim(rotMax,1);
clear rotMax;

rotMean = mean(rot);
fm(:,:,12) = shiftdim(rotMean,1);
clear rotMean;

rotVar = var(rot);
fm(:,:,13) = shiftdim(rotVar,1);
clear rotVar;

rotMedian = median(rot);
fm(:,:,14) = shiftdim(rotMedian,1);
clear rotMedian;

%rotMax - rotMin
fm(:,:,15) = fm(:,:,11) - fm(:,:,10);

clear rot;

%% Histogram
% window size is similar to orientatinFilter given above
csHalf = floor(csHist/2);
tic;
for i=1:length(imIn(:))
[r,c] = ind2sub(size(imIn),i);
if r < oriFiltLen | c < oriFiltLen | r > size(imIn,1)-oriFiltLen | c > size(imIn,2)-oriFiltLen
  continue
else
sub = imIn(r-csHalf:r+csHalf, c-csHalf:c+csHalf);
sub = norm01(sub) * 100;
[m,v,h] = meanvar(sub);
fm(r,c,16:25) = h;
fm(r,c,26) = m;
fm(r,c,27) = v;
end
end
toc;

[eig1, eig2, cw] = structureTensorImage2(imIn, 1, 1);
[gx,gy,mag] = gradientImg(imIn,1);
clear gx
clear gy

%% 60 more features - DOG, gradients, smoothened
tic;
%used to be 20
k = 0;
for i=1:1:10
ims = imsmooth(double(imIn),i);
k = k + 1;
fm(:,:,27+k) = ims;
k = k + 1;
fm(:,:,27+k) = imsmooth(eig1./eig2,i);
k = k + 1;
fm(:,:,27+k) = imsmooth(mag,i);
[gx,gy,mags] = gradientImg(imIn,i);
k = k + 1;
fm(:,:,27+k) = mags;

    for j=1:2:i-2
        k = k + 1;
        fm(:,:,27+k) = ims - single(imsmooth(double(imIn),j));
    end

end
toc;
% feature count = 60 + 27 = 87
%% DOG
smoothFM = fm(:,:,78:87);
fm(:,:,88) = shiftdim(var(shiftdim(smoothFM,2)),1);
fm(:,:,89) = norm01(imsmooth(double(imIn),2) - ...
         imsmooth(double(imIn),50));
% fm(:,:,90) = imIn;

disp('getFeatures_nonMembrane.m: Feature extraction finished for single image.')