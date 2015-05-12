function cocMat = getXcorrZYstackZ(inputImageStackFileName,maxShift,minShift,maxNumImages)
% calculate the c.o.c of ZY face along the Z axis

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
[numR,numImages,numC] = size(inputImageStack); 

% TODO: current we take the first n images for the estimation. Perhaps we
% can think of geting a random n images.
disp('Estimating similarity curve using correlation coefficient of shifted XY sections ...')
if(maxNumImages>numImages)
    maxNumImages = numImages;
    disp('maxNumImages > numImages. using numImages = %d instead',numImages);
end
numShifts = maxShift - minShift + 1;
cocMat = zeros(maxNumImages,numShifts);
for z=1:maxNumImages 
    k=0;
    for g=minShift:maxShift
        A = zeros(numR,numC-g);
        B = zeros(numR,numC-g);

        A(:,:) = inputImageStack(:,z,1+g:numC);
        B(:,:) = inputImageStack(:,z,1:numC-g);
        k=k+1;
        cocMat(z,k) = corr2(A,B);
    end
end

%% plot
% titleStr = 'Coefficient of Correlation using XY sections along X axis';
% xlabelStr = 'Shifted pixels';
% ylabelStr = 'Coefficient of Correlation';
% transparent = 0;
% shadedErrorBar((1:maxShift),mean(cocMat,1),std(cocMat),'g',transparent,...
%     titleStr,xlabelStr,ylabelStr);
