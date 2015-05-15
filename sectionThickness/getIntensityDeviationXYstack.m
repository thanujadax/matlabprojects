function sigmaMat = getIntensityDeviationXYstack...
    (inputImageStackFileName,maxShift,minShift,maxNumImages)
% calculate the sd of intensity difference along the xy plane. i.e. parallel to the
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
if(maxNumImages>numImages)
    maxNumImages = numImages;
    str1 = sprintf('maxNumImages > numImages. using numImages = %d instead',numImages);
    disp(str1)
end
numShifts = maxShift - minShift + 1;
sigmaMat = zeros(maxNumImages,numShifts);

% TODO: current we take the first n images for the estimation. Perhaps we
% can think of geting a random n images.
disp('Estimating similarity curve using SD of intensity differences across shifted XY sections')

for z=1:maxNumImages
    I = inputImageStack(:,:,z);
    [numR,numC] = size(I);
    k = 0;
    for g=minShift:maxShift
        d1I = (I(1+g:size(I,1),:)-I(1:size(I,1)-g,:));
        d2I = (I(:,1+g:size(I,2))-I(:,1:size(I,2)-g));
        k = k + 1;
        sigmaMat(z,k) = std([d1I(:);d2I(:)]);
    end
end

%% plot
% titleStr = 'SD of pixel intensity deviations using shifted XY sections';
% xlabelStr = 'Shifted pixels';
% ylabelStr = 'SD of pixel intensity deviation';
% transparent = 1;
% shadedErrorBar((1:maxShift),mean(sigmaMat,1),std(sigmaMat),'g',transparent,...
%     titleStr,xlabelStr,ylabelStr);