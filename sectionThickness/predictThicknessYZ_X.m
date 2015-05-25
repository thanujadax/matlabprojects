function predictedThickness = predictThicknessYZ_X...
        (inputImageStackFileName,meanVector,inputResolution,...
        distMin,method)
    
% read image stack
% calculate pairwise c.o.c of each adjacent pair of images YZ_X
% interpolate the decay curve to predict thickness

inputImageStack = readTiffStackToArray(inputImageStackFileName);

[sizeR,sizeC,sizeZ] = size(inputImageStack);
distMax = numel(meanVector);
A = zeros(sizeR,sizeZ);
B = zeros(sizeR,sizeZ);

predictedThickness = zeros((sizeC-1),1);

for i=1:sizeC-1
    A(:,:) = inputImageStack(:,i,:);
    B(:,:) = inputImageStack(:,(i+1),:);
    coc = corr2(A,B);
    predictedThickness(i) = interp1(meanVector,(distMin:distMax-1),coc,method)...
                            .* inputResolution;
end