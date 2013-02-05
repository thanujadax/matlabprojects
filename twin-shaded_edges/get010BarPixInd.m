function [pixelIndUp pixIndDown] = get010BarPixInd(r,c,orientation,barLength,barWidth,offWidth,numRows,numCols)
% works for any orientation given in degrees
% returns the indices of the pixels in the top part of the bar which
% includes pixel (r,c)

% define bar with 1,1 as one corner
barOn = ones(barWidth,barLength);
% get the center of the bar. 
c0 = floor((barLength+1)/2);
r0 = floor((barWidth+1)/2) + offWidth;

% % set the lower half of the bar to -1
% lowerRowStart = r0+1;
% lowerRowEnd = barWidth;
% barOn(lowerRowStart:lowerRowEnd,:) = -1;

% create 010 kernel
barOff = ones(offWidth,barLength) .* (-1);
barOn = [barOff;barOn;barOff];

% rotate bar around it's mid point
if(orientation>0)
    orientation = 180 - orientation;
    barOn = imrotate(barOn,orientation);
end

[rowsUp colsUp] = find(barOn>0); % positive pixels (middle)
[rowsDown colsDown] = find(barOn<0);
% shift the bar to be centered around (r,c)
rShift = r - r0;
cShift = c - c0;

% (+1)
rowsUp = rowsUp + rShift;
colsUp = colsUp + cShift;
pixelIndUp = sub2ind([numRows numCols],rowsUp,colsUp);

% bottom half (-1)
rowsDown = rowsDown + rShift;
colsDown = colsDown + cShift;
pixIndDown = sub2ind([numRows numCols],rowsDown,colsDown);