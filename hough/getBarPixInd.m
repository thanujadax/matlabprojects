function pixelInd = getBarPixInd(r,c,orientation,barLength,barWidth,numRows,numCols)
% works for any orientation given in degrees
% gives the indices of the pixels that form a bar of given length and with
% and orientation, centered around the given point (r,c)

% define bar with 1,1 as one corner
barOn = ones(barWidth,barLength);
% get the center of the bar
c0 = floor((barLength+1)/2);
r0 = floor((barWidth+1)/2);
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

% check if there are any elements of 'rows' which is greater than the value
% of 'numRows'. avoid this situation since it will give a runtime error.
% the same for 'cols'.
if(numel(find(rows>numRows))>0 || numel(find(cols>numCols))>0)
    pixelInd = -1;
else
    pixelInd = sub2ind([numRows numCols],rows,cols);
end
% if(numel(find(rows>numRows))==0 || numel(find(cols>numCols))==0 || ...
%         find(rows<1))
%     pixelInd = sub2ind([numRows numCols],rows,cols);
% else
%     pixelInd = -1;
% end