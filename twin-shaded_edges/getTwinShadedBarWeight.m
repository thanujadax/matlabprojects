function pixWeight = getTwinShadedBarWeight(img,r,c,orientation,barLength,barWidth)
% returns the cumulative weight of all the pixels inside the oriented bar
% centered at (r,c) in img

% the bar is divided into two parts. The top is bright and contains the
% pixel in focus (r,c). The bottom is dark. This means when we count the
% support of the pixels, the top part counts the support normally and then
% support of the bottom part is subtracted from it.

[numRows numCols] = size(img);

[pixIndUp pixIndDown] = getTwinBarPixInd(r,c,orientation,barLength,barWidth,numRows,numCols);
% pixIndUp contain the indices of the pixels that are bright (membrane)
% pixIndDown contain the indices of the pixels that are dark (neuron)

pixValsUp = img(pixIndUp); % should contain mostly positive values at the edges
pixValsDown = img(pixIndDown); % should contain mostly negative values at the edges

pixWeight = sum(pixValsUp) - sum(pixValsDown);
