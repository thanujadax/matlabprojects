function [myelinSegments] = getMyelinSegments1(...
                    probabilityMap,params)

% TODO: description

% Inputs
%     probabilityMap - myelin pixel probabilities  
%     params.sigmaGaussianBlur = 2;
%     params.maskSizeGaussianBlur = 7;
%     params.threshold_pixelIntensity = 0.15; % to generate binary images for myelin labels
%     params.numPixelsCC = 3000; % minimum number of pixels for CC to be considered as myelin

% Output
%   myelinSegments - binary image where 1 indicates myelin and 0 otherwise

% gaussian blur
gbImage = gaussianFilter(probabilityMap,params.sigmaGaussianBlur,...
                params.maskSizeGaussianBlur);
% binary thresholding
thresholdedGbImage = im2bw(gbImage,params.threshold_pixelIntensity);

% extraction of connected components
CC = bwconncomp(thresholdedGbImage);
numPixels = cellfun(@numel,CC.PixelIdxList);
ccInd = find(numPixels>params.numPixelsCC);

myelinSegments = zeros(size(probabilityMap));
if(~isempty(ccInd))
    for i=1:length(ccInd)
        pixIndLogical_i= CC.PixelIdxList{ccInd(i)};
        myelinSegments(pixIndLogical_i) = 1;
    end
end