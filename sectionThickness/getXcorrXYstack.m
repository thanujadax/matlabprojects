function xcorrMat = getXcorrXYstack(inputImageStackFileName,maxShift,minShift,maxNumImages)
% calculate the correlation of XY face along the Y axis . i.e. parallel to the
% cutting plane where we have maximum resolution (5nmx5nm for FIBSEM)

% Inputs:
% imageStack - image stack (tif) for which the thickness has to be
% estimated. This has to be registered along the z axis already.

inputImageStack = readTiffStackToArray(inputImageStackFileName);
% inputImageStack is a 3D array where the 3rd dimension is along the z axis

% estimate the correlation curve (mean and sd) from different sets of
% images within the max window given by maxShift

% initially use just one image

% I = double(imread(imageStack));
numImages = size(inputImageStack,3);

% TODO: current we take the first n images for the estimation. Perhaps we
% can think of geting a random n images.
disp('Estimating similarity curve using correlation coefficient of shifted XY sections ...')
if(maxNumImages>numImages)
    maxNumImages = numImages;
    disp('maxNumImages > numImages. using numImages = %d instead',numImages);
end
numShifts = maxShift - minShift + 1;
xcorrMat = zeros(maxNumImages,numShifts);

for z=1:maxNumImages
    I = inputImageStack(:,:,z);
    [numR,numC] = size(I);
    k = 0;
    for g=minShift:maxShift
        k = k + 1;
        A = zeros(numR-g,numC);
        B = zeros(numR-g,numC);

        A(:,:) = I(1+g:size(I,1),:);
        B(:,:) = I(1:size(I,1)-g,:);
        
        xcorrMat(z,k) = corr2(A,B);
    end
end

%% plot
% titleStr = 'Coefficient of Correlation using XY sections along Y axis';
% xlabelStr = 'Shifted pixels';
% ylabelStr = 'Coefficient of Correlation';
% transparent = 0;
% shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g',transparent,...
%     titleStr,xlabelStr,ylabelStr);
