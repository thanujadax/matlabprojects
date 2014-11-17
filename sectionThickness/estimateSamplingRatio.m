function [g,m12] = estimateSamplingRatio(I,maxG)
% ESTIMATESAMPLINGRATIO estimates distance.
%
% from Sporring 2014

% I - 2 images packed as NxMx2
% maxG - maximum sampling distance to evaluate

normalize = @(I) (I-mean(I(:)))/std(I(:));
I(:,:,1) = normalize(I(:,:,1));
I(:,:,2) = normalize(I(:,:,2));
m12 = zeros(1,maxG);

for g = 1:maxG
    d1I = (I(1+g:size(I,1),:,:)-I(1:size(I,1)-g,:,:));
    d2I = (I(:,1+g:size(I,2),:)-I(:,1:size(I,2)-g,:));
    m12(g) = std([d1I(:);d2I(:)]);
end

d3I = I(:,:,2)-I(:,:,1);
m3 = std(d3I(:));
g = interp1(m12,1:maxG,m3);
