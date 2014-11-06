function overlayImage = overlayWatershedOnImage(originalImg,watershedImg)
% overlay watershed edges on to original image


watershedEdgePixels = (watershedImg==0);
figure();
watershedColoredEdges = cat(3,watershedEdgePixels,...
    zeros(size(watershedImg)),zeros(size(watershedImg))); % red
imshow(watershedColoredEdges);
hold on
h = imshow(originalImg);
hold off

set(h, 'AlphaData', watershedEdgePixels);
