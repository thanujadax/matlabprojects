function xcorrMat = getXcorrShiftedImg(imageFileName,maxShift)
% shifts the image along the y axis and calculates the correlation
% coefficient

I = double(imread(imageFileName));
[numR,numC] = size(I);
xcorrMat = zeros(1,maxShift);

for g=1:maxShift
    A = zeros(numR-g,numC);
    B = zeros(numR-g,numC);
    
    A(:,:) = I(1+g:size(I,1),:);
    B(:,:) = I(1:size(I,1)-g,:);
    
    xcorrMat(g) = corr2(A,B);
end