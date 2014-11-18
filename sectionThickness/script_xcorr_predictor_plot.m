% script to plot estimator(s) of section thickness with their confidence
% (or variance)

% estimator based on correlation coefficient
xcorrMeasure = 1;
% if 0, use sd of intensity difference method

maxShift = 15;



% read set of input images
inputImageDir = '/home/thanuja/projects/inputData/trainingHalf/raw';
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