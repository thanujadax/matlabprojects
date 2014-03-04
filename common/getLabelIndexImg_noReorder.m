function [labelIndImg0,numLabels] = getLabelIndexImg_noReorder(labelImage)

% labelImage is RGB (mxnx3 matrix). 
% labelIndImg0 replaces colors with unique numbers so that the image
% becomes a 2D matrix

grayImg = rgb2gray(labelImage);
labelIndImg0 = im2uint16(grayImg);

% re-index from 1 to n, in order
labels_old = unique(labelIndImg0);
labels_old = labels_old(labels_old>0);

numLabels = numel(labels_old);
