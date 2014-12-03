function [t_median,t_mean,t_var] = getThicknessForRegisteredPairs...
    (image1_patches,image2_patches,thicknessCurve,maxThicknessPix)

% given registered patches from two images, estimate the thickness of
% (distance between) the original image pair

% Input:
%   thicknessCurve - pixel gap vs correlation curve

% Outptut:
    

image1_patches = normalize(image1_patches);
image2_patches = normalize(image2_patches);

[~,~,numPatches] = size(image1_patches);
thickness = zeros(numPatches,1);

for i=1:numPatches
    c = corr2(image1_patches(:,:,i),image2_patches(:,:,i));
    thickness(i) = interp1(thicknessCurve,1:maxThicknessPix,c);
end

t_median = median(thickness);
t_mean = mean(thickness);
t_var = var(thickness);