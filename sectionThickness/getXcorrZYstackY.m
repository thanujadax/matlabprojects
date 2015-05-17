function xcorrMat = getXcorrZYstackY(inputImageStackFileName,maxShift,minShift,maxNumImages)
% calculate c.o.c curve on the ZY plane, shifting the (same) image along
% the Y axis. Multiple images are used to get an average estimate of the
% decay of c.o.c.

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

% TODO: current we take the first n images for the estimation. Perhaps we
% can think of geting a random n images.
disp('Estimating similarity curve using zy sections, shifting along Y ...')
if(maxNumImages>numImages)
    maxNumImages = numImages;
    str1 = sprintf('maxNumImages > numImages. using numImages = %d instead',numImages);
    disp(str1)
end
numShifts = maxShift - minShift + 1;
xcorrMat = zeros(maxNumImages,numShifts);
for z=1:maxNumImages
    k=0;
    for g=minShift:maxShift
        A = zeros(numR-g,numC);
        B = zeros(numR-g,numC);   
        A(:,:) = inputImageStack(1+g:size(inputImageStack,1),z,:);
        B(:,:) = inputImageStack(1:size(inputImageStack,1)-g,z,:);  % with shift
        k=k+1;
        xcorrMat(z,k) = corr2(A,B);
    end
end
%% plot
% titleStr = 'Coefficient of Correlation using ZY sections, along Y';
% xlabelStr = 'Shifted pixels';
% ylabelStr = 'Coefficient of Correlation';
% transparent = 0;
% shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g',transparent,...
%     titleStr,xlabelStr,ylabelStr);