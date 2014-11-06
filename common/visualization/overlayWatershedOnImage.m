function watershedColoredEdges = overlayWatershedOnImage(originalImg,watershedImg)
% overlay watershed edges on to original image


% watershedEdgePixels = (watershedImg==0);
% watershedColoredEdges = cat(3,watershedEdgePixels,...
%     zeros(size(watershedImg)),zeros(size(watershedImg))); % red
% 
% imshow(originalImg);
% hold on
% h = imshow(watershedColoredEdges);
% set(h, 'AlphaData', 0.5);
% hold off




watershedEdgePixels = (watershedImg==0);
watershedColoredEdges = cat(3,watershedEdgePixels,...
    zeros(size(watershedImg)),zeros(size(watershedImg))); % red
% figure;imshow(originalImg)
% figure;imshow(watershedColoredEdges);
figure;imshowpair(originalImg,watershedColoredEdges)