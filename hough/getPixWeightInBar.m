function pixWeight = getPixWeightInBar(img,r,c,orientation,barLength,barWidth)
% returns the cumulative weight of all the pixels inside the oriented bar
% centered at (r,c) in img

[numRows numCols] = size(img);
pixInd = getBarPixInd(r,c,orientation,barLength,barWidth,numRows,numCols);
pixVals = img(pixInd);
pixWeight = sum(pixVals);