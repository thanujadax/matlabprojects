% script to plot estimator(s) of section thickness with their confidence
% (or variance)

xyResolution = 5; % nm
% estimator based on correlation coefficient
xcorrMeasure = 0;
% if 0, use sd of intensity difference method

maxShift = 15;



% read set of input images
inputImageDir = '/home/thanuja/projects/inputData/trainingHalf/raw';
% inputImageDir = '/home/thanuja/Dropbox/data/fibsem/smallSet2/raw';
allImageFiles = dir(fullfile(inputImageDir,'*.png'));

% calculate cross correlation for different pixel intervals
% for each file
numImg = length(allImageFiles);

xcorrMat = zeros(numImg,maxShift);

for i=1:numImg
    disp(i);
    imageFileName = fullfile(inputImageDir,allImageFiles(i).name);
    if(xcorrMeasure)
        xcorrMat(i,:) = getXcorrShiftedImg(imageFileName,maxShift);
    else
        xcorrMat(i,:) = getIntensityDeviationShiftedImg(imageFileName,maxShift);
    end
end

% plot
shadedErrorBar((1:maxShift),mean(xcorrMat,1),std(xcorrMat),'g');

% predict section thickness for the data set
relZresolution = zeros(1,numImg-1); % relative to xy pix resolution
for i = 1:numImg-1
   imageName1 = fullfile(inputImageDir,allImageFiles(i).name);
   imageName2 = fullfile(inputImageDir,allImageFiles(i+1).name);
   if(xcorrMeasure)
       relZresolution(i) = getRelativeDistance_cc(imageName1,imageName2,mean(xcorrMat,1),maxShift);
   else
       relZresolution(i) = getRelativeDistance_intensitydiffsd(imageName1,imageName2,mean(xcorrMat,1),maxShift);
   end
end

estimatedSectionThickness = relZresolution .* xyResolution;