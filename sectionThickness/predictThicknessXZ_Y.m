function predictedThickness = predictThicknessXZ_Y...
        (inputImageStackFileName,meanVector,inputResolution,...
        distMin,method)
    
% read image stack
% calculate pairwise c.o.c of each adjacent pair of images YZ_X
% interpolate the decay curve to predict thickness

inputImageStack = readTiffStackToArray(inputImageStackFileName);

[sizeR,sizeC,sizeZ] = size(inputImageStack);
distMax = numel(meanVector);
A = zeros(sizeC,sizeZ);
B = zeros(sizeC,sizeZ);

predictedThickness = zeros((sizeR-1),1);

for i=1:sizeR-1
    A(:,:) = inputImageStack(i,:,:);
    B(:,:) = inputImageStack((i+1),:,:);
    coc = corr2(A,B);
    predictedThickness(i) = interp1(meanVector,(distMin:distMax-1),coc,method)...
                            .* inputResolution;
end