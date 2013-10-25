function [labelIndImg,numLabels] = getLabelIndexImg(labelImage)

grayImg = rgb2gray(labelImage);
labelIndImg0 = im2uint16(grayImg);

% re-index from 1 to n, in order
labels_old = unique(labelIndImg0);
labels_old = labels_old(labels_old>0);

sortedLabels = sort(labels_old);

numLabels = numel(sortedLabels);

labelIndImg = labelIndImg0;

for i=1:numLabels
    pixelsForLabel_i = (labelIndImg0==sortedLabels(i));
    labelIndImg(pixelsForLabel_i) = i;
end
