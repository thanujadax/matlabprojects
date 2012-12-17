function patchLines = localHoughLines(localHoughSpaces,R,T,imgIn,bb,maxLinesPerPatch,...
            peakThresh,houghSupNHood,fillGap,minLength,maxHoughPeak)
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
% maxHoughPeak
        
        
% output:
% patchLines - a cell array containing the lines of each patch
%   

% get the size of localHoughSpaces cell array
[rows cols] = size(localHoughSpaces);
[imgRows imgCols] = size(imgIn);

% initialize the struct array patchLines
patchLines = cell(rows,cols);
 %lineStruct = struct('point1',[],'point2',[],'theta',nan,'rho',nan);
 %patchLines(rows,cols) = lineStruct;
        

for i = 1:rows
    startRow = (i-1)*bb +1;
    if(startRow==imgRows)
        break;
    end
    stopRow = startRow + bb -1;
    if(stopRow>imgRows)
        stopRow = imgRows;
    end
    for j = 1:cols
        % for the H of the (i,j)th patch 
        P  = houghpeaks(localHoughSpaces{i,j},maxLinesPerPatch,'threshold',ceil(peakThresh*maxHoughPeak),...
            'NHoodSize',houghSupNHood);
        
        startCol = (j-1)*bb +1;       
        stopCol = startCol + bb -1;
        if(startCol==imgCols)
            break
        end
        if(stopCol<=imgCols)
            patchLines{i,j} = houghlines(imgIn(startRow:stopRow, startCol:stopCol)...
                ,T,R,P,'FillGap',fillGap,'MinLength',minLength);
        else
            stopCol = imgCols;
            patchLines{i,j} = houghlines(imgIn(startRow:stopRow, startCol:stopCol)...
                ,T,R,P,'FillGap',fillGap,'MinLength',minLength);
        end
        
    end

end