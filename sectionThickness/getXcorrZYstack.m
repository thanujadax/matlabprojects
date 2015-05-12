function xcorrMat = getXcorrZYstack(inputImageStackFileName,maxShift,maxNumImages)
% calculate the correlation of the ZY plane along the X axis.

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
numImages = size(inputImageStack,2); % x axis

A = zeros(numR,numC);
B = zeros(numR,numC);

z = 1; % starting image
% TODO: current we take the first n images for the estimation. Perhaps we
% can think of geting a random n images.
disp('Estimating similarity curve using zy sections ...')
numImages = numImages - maxShift;
if(maxNumImages>numImages)
    maxNumImages = numImages;
    disp('maxNumImages > numImages. using numImages = %d instead',numImages);
end
xcorrMat = zeros(maxNumImages,maxShift);

for z=1:maxNumImages
    for g=1:maxShift
        A(:,:) = inputImageStack(:,z,:);
        B(:,:) = inputImageStack(:,z+g,:);  % with shift
        xcorrMat(z,g) = corr2(A,B);
    end
end
%% plot
titleStr = 'Coefficient of Correlation using ZY sections';
xlabelStr = 'Shifted pixels';
ylabelStr = 'Coefficient of Correlation';
transparent = 0;
shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g',transparent,...
    titleStr,xlabelStr,ylabelStr);