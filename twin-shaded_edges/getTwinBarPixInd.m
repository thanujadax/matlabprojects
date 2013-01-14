function [pixelIndUp pixIndDown] = getTwinBarPixInd(r,c,orientation,barLength,barWidth,numRows,numCols)
% works for any orientation given in degrees
% returns the indices of the pixels in the top part of the bar which
% includes pixel (r,c)

% define bar with 1,1 as one corner
barOn = ones(barWidth,barLength);
% get the center of the bar. Here, since the width is even, the center would
% be biased to one side. i.e. the top part.
c0 = floor((barLength+1)/2);
r0 = floor((barWidth+1)/2);

% set the lower half of the bar to -1
lowerRowStart = r0+1;
lowerRowEnd = barWidth;
barOn(lowerRowStart:lowerRowEnd,:) = -1;

% rotate bar around it's mid point
if(orientation>0)
    orientation = 180 - orientation;
    barOn = imrotate(barOn,orientation);
end

[rowsUp colsUp] = find(barOn); % positive pixels (up)
[rowsDown colsDown] = find(barOn<0);
% shift the bar to be centered around (r,c)
rShift = r - r0;
cShift = c - c0;

% top half (+1)
rowsUp = rowsUp + rShift;
colsUp = colsUp + cShift;
pixelIndUp = sub2ind([numRows numCols],rowsUp,colsUp);

% bottom half (-1)
rowsDown = rowsDown + rShift;
colsDown = colsDown + cShift;
pixIndDown = sub2ind([numRows numCols],rowsDown,colsDown);
