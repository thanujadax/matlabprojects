function relZresolution = predictThicknessFromCurveFromMultiplePairs(...
        imageStackFileName,xcorrMat,minShift,maxShift,calibrationMethod,numPairs)
% Returns the section thickness relative to the xy resolution. Multiply by
% xyResolution to get the actual thickness.

% calibrationMethod
% % 1 - c.o.c across XY sections, along X
% % 2 - c.o.c across XY sections, along Y axis
% % 3 - c.o.c across ZY sections, along x axis
% % 4 - c.o.c across ZY sections along Y
% % 5 - c.o.c acroxx XZ sections, along X
% % 6 - c.o.c acroxx XZ sections, along Y
% % 7 - c.o.c across XY sections, along Z
% % 8 - c.o.c across ZY sections, along Z
% % 9 - c.o.c. across XZ sections, along Z
% % 10 - SD of XY per pixel intensity difference
% TODO: methods robust against registration problems

% numPairs -
%   if set to 1, use sections i and i + 1 to estimate thickness of i
%   if 2, use sections i, i+1 and i+2 to estiamte the thickness of i
    
inputImageStack = readTiffStackToArray(imageStackFileName);
numImg = size(inputImageStack,3);

relZresolution = zeros(numPairs,numImg-1); % relative to xy pix resolution
str1 = sprintf('Calculating distances using %d pairs',numPairs);
disp(str1);

if(numPairs > 0)
    if(calibrationMethod==10)
        for i = 1:numImg-1
            image1 = inputImageStack(:,:,i);
            image2 = inputImageStack(:,:,(i+1));
            % calculate the distance between two images based on the SD of
            % pixel differences
            deviationSigma = getPixIntensityDeviationSigma(image1,image2);
            relZresolution(1,i) = getRelativeDistance(mean(xcorrMat,1),minShift,maxShift,deviationSigma);
        end
    else
        for i = 1:numImg-1
            image1 = inputImageStack(:,:,i);
            image2 = inputImageStack(:,:,(i+1));
            % calculate the distance between the two images based on the
            % correlation coefficient
            relZresolution(1,i) = getRelativeDistance_cc2(image1,image2,mean(xcorrMat,1),maxShift,minShift);
        end
    end
end

if(numPairs>1)
    if(calibrationMethod==10)
        for i = 1:numImg-2
            image1 = inputImageStack(:,:,i);
            image2 = inputImageStack(:,:,(i+2));
            % calculate the distance between two images based on the SD of
            % pixel differences
            deviationSigma = getPixIntensityDeviationSigma(image1,image2);
            relZresolution(2,i) = getRelativeDistance(mean(xcorrMat,1),maxShift,minShift,deviationSigma);
        end
    else
        for i = 1:numImg-2
            image1 = inputImageStack(:,:,i);
            image2 = inputImageStack(:,:,(i+2));
            % calculate the distance between the two images based on the
            % correlation coefficient
            relZresolution(2,i) = getRelativeDistance_cc2(image1,image2,mean(xcorrMat,1),maxShift,minShift);
        end
    end    
end