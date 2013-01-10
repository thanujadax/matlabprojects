function pixelInd = getBarPixInd(r,c,orientation,barLength,barWidth,numRows,numCols)

% define bar with 1,1 as one corner
barOn = ones(barWidth,barLength);
% get the center of the bar
c0 = floor((barLength-1)/2)+1;
r0 = floor((barWidth-1)/2)+1;
% rotate bar around it's mid point
if(orientation>0)
    orientation = 180 - orientation;
    barOn = imrotate(barOn,orientation);
end
% 'crop' makes sure that the rotated bar is the same size as 
[rows cols] = find(barOn);
% shift the bar to be centered around (r,c)
rShift = r - r0;
cShift = c - c0;
%[rows cols] = ind2sub([barWidth barLength],barInd);
rows = rows + rShift;
cols = cols + cShift;
pixelInd = sub2ind([numRows numCols],rows,cols);