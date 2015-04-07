function xcorrMat = getXcorrZYstack(inputImageStackFileName,maxShift,maxNumImages)
% calculate the correlation along the zy plane. i.e. perpendicular to the
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
numR = size(inputImageStack,1);
numC = size(inputImageStack,3); % z axis

xcorrMat = zeros(maxNumImages,maxShift);

A = zeros(numR,numC);
B = zeros(numR,numC);
z = 1;

for z=1:maxNumImages
    for g=1:maxShift
        A(:,:) = InputImageStack(:,z,:);
        B(:,:) = InputImageStack(:,z+g,:);  % with shift
        xcorrMat(z,g) = corr2(A,B);
    end
end