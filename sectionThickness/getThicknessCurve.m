function [thicknessCurve, t_std] = getThicknessCurve...
                                (inputImageDir,imageType,maxShift)

% Inputs:
%   inputImageDir
%   maxShift - max shift in pixels to consider for the curve

% Outputs:
%   thicknessCurve - correlation coefficient vs pix gap  


allImageFiles = dir(fullfile(inputImageDir,strcat('*.',imageType)));
numImg = length(allImageFiles);

xcorrMat = zeros(numImg,maxShift);

for i=1:numImg
    imageFileName = fullfile(inputImageDir,allImageFiles(i).name);
    if(xcorrMeasure)
        xcorrMat(i,:) = getXcorrShiftedImg(imageFileName,maxShift);
    else
        xcorrMat(i,:) = getIntensityDeviationShiftedImg(imageFileName,maxShift);
    end
end

thicknessCurve = mean(xcorrMat,2);
t_std = std(xcorrMat,0,2);