function patchLines = localHoughLines(localHoughSpaces,R,T,imgIn,bb,maxLinesPerPatch,...
            peakThresh,houghSupNHood,fillGap,minLength)
% inputs:
% localHoughSpaces -
% R - rho values
% T - theta values
% imgIn - original image
% bb - patch size
% maxLinesPerPatch -  
% peakThresh - fraction of max(H) to be used for peak selection
% houghSupNHood - 
% fillGap - if there's a gap less than this between 2 line segments, merge
% them
% minLength - minimum line length to consider
        
        
% output:
% patchLines - a cell array containing the lines of each
% patch. Dimenstions : [numRowPatches numColPatches]
%   patchLines.start = [r1 c1]
%   patchLines.ends = [r2 c2]

% get the size of localHoughSpaces cell array
[rows cols] = size(localHoughSpaces);

% initialize the cell array patchLines
patchLines = cell(rows,cols);

for i = 1:rows
    for j = 1:cols
        % for the H of the (i,j)th patch 
        P  = houghpeaks(localHoughSpaces{i,j},maxLinesPerPatch,'threshold',ceil(peakThresh*max(H(:))),...
            'NHoodSize',houghSupNHood);
        startRow = (i-1)*bb +1;
        startCol = (j-1)*bb +1;
        stopRow = startRow + bb -1;
        stopCol = startCol + bb -1;
        patchLines{i,j} = houghlines(imgIn(startRow:stopRow, startCol:stopCol)...
            ,T,R,P,'FillGap',fillGap,'MinLength',minLength);
        
    end
end