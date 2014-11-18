function [g,m12] = estimateRelResolution_xcorr(I,maxG)
% based on correlation coefficient

% I - 2 images packed as NxMx2
% maxG - maximum sampling distance to evaluate

normalize = @(I) (I-mean(I(:)))/std(I(:));
I(:,:,1) = normalize(I(:,:,1));
I(:,:,2) = normalize(I(:,:,2));
m12 = zeros(1,maxG);

[numR,numC,~] = size(I);

for g = 1:maxG
    
    A = zeros(numR-g,numC);
    B = zeros(numR-g,numC);
    
    A(:,:) = I(1+g:size(I,1),:,2);
    B(:,:) = I(1:size(I,1)-g,:,2);
    
    m12(g) = corr2(A,B);
    
    
    % d2I = (I(:,1+g:size(I,2),:)-I(:,1:size(I,2)-g,:));
    % m12(g) = std([d1I(:);d2I(:)]);
    % autocorr_xy = corr2(d1I,d2I);
end

% d3I = I(:,:,2)-I(:,:,1);
% m3 = std(d3I(:));
m3 = corr2(I(:,:,1),I(:,:,2));
g = interp1(m12,1:maxG,m3);