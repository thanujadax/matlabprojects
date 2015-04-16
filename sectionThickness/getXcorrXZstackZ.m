function cocMat = getXcorrXZstackZ(inputImageStackFileName,maxShift,maxNumImages)
% calculate the correlation of XZ face along the Z axis.

% Inputs:
% imageStack - image stack (tif) for which the thickness has to be
% estimated. This has to be registered along the z axis already.

% Outputs:
%   cocMat - Matrix containg coefficient of correlation for image pairs.
%   Each row corresponds to a different starting image. The distance
%   increases with the column index

inputImageStack = readTiffStackToArray(inputImageStackFileName);
% inputImageStack is a 3D array where the 3rd dimension is along the z axis

% estimate the correlation curve (mean and sd) from different sets of
% images within the max window given by maxShift

% initially use just one image

% I = double(imread(imageStack));
numR = size(inputImageStack,1);
numC = size(inputImageStack,3); % z axis

cocMat = zeros(maxNumImages,maxShift);

% TODO: current we take the first n images for the estimation. Perhaps we
% can think of geting a random n images.
disp('Estimating similarity curve using correlation coefficient of shifted XY sections ...')
for z=1:maxNumImages
    I = inputImageStack(z,:,:);
    [numR,numC] = size(I);
    for g=1:maxShift
        A = zeros(numR-g,numC);
        B = zeros(numR-g,numC);

        A(:,:) = I(1+g:size(I,1),:);
        B(:,:) = I(1:size(I,1)-g,:);
        
        cocMat(z,g) = corr2(A,B);
    end
end

%% plot
titleStr = 'Coefficient of Correlation using XZ sections along Z axis';
xlabelStr = 'Shifted pixels';
ylabelStr = 'Coefficient of Correlation';
transparent = 0;
shadedErrorBar((1:maxShift),mean(cocMat,1),std(cocMat),'g',transparent,...
    titleStr,xlabelStr,ylabelStr);