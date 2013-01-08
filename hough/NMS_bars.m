function voteMat = NMS_bars(voteMat,orientation,barLength,barWidth)

% non-maximal suppression
% for each non zero pixel, make it's value zero if there's at least one
% value equal or more than it. 

% output:
%   voteMat - the vote matrix for the input orientation after NMS

% inputs:
%   voteMat - vote matrix for the current orientation after thresholding
%   orientation
%   barLength
%   barWidth

maxVote = max(max(voteMat));

[r,c,vote] = find(voteMat);

numPixels = numel(r);

for j=1:numPixels
    % get neighbor hood of (r(i),c(i))
    if(orientation==0)
        % horizontal bar
        startRow = r - ceil(barWidth/2);
        stopRow = startRow + barWidth;
        startCol = c - ceil(barLength/2);
        stopCol = startCol + barLength;

        [pxlR, pxlC] = find(voteMat>=vote(j)); % pixels having a higher vote    

        rowIndInd = intersect(find(pxlR>startRow),find(pxlR<stopRow));    
        colIndInd = intersect(find(pxlC>startCol),find(pxlC<stopCol));    

        pixInd = intersect(rowIndInd,colIndInd);
        
        if(numel(pixInd)>0)
            % this point is not a maximum. set it to zero
            voteMat(r(j),c(j)) = 0;
        end
        
    elseif(orientation==90)
        % vertical bar
        startRow = r - ceil(barLength/2);
        stopRow = startRow + barLength;
        startCol = c - ceil(barWidth/2);
        stopCol = startCol + barWidth;

        [pxlR, pxlC] = find(voteMat>=vote(j)); % pixels having a higher vote   

        rowIndInd = intersect(find(pxlR>startRow),find(pxlR<stopRow));    
        colIndInd = intersect(find(pxlC>startCol),find(pxlC<stopCol));  

        pixInd = intersect(rowIndInd,colIndInd);   
        
        if(numel(pixInd)>0)
            % this point is not a maximum. set it to zero
            voteMat(r(j),c(j)) = 0;
        end   
        
    elseif(orientation==45)
        % get the list of all pixels in the bar
        % 1: the line through (c,r)    
        b = r-c;        % b is the y-intercept of the straight line y=ax+b going through (c,r)
        cStart = c - floor(cosd(orientation)*barLength/2);
        cEnd = c + floor(cosd(orientation)*barLength/2);

        cRange = cStart:cEnd;
        rRange = cRange + b;    % from the eqn of the straight line

        line0 = [rRange' cRange'];
        % for each additional line
        nLinePairs = (barWidth-1)/2; % number of additional pairs of lines
        for i = 1:nLinePairs
            rUp = r-i;
            bUp = rUp-c;
            cStartUp = cStart+i;
            cEndUp = cEnd+i;
            cRangeUp = cStartUp:cEndUp;
            rRangeUp = cRangeUp + bUp;

            rDown = r+i;
            bDown = rDown-c;
            cStartDown = cStart-i;

            cEndDown = cEnd-i;
            cRangeDown = cStartDown:cEndDown;
            rRangeDown = cRangeDown + bDown;

            lineUp = [rRangeUp' cRangeUp'];
            lineDown = [rRangeDown' cRangeDown'];

            line0 = [line0; lineUp; lineDown];        
        end
        % line0 contains (r,c) of each point in the bar as a column vector
        % count the pixels in line0 above the vote of current pixel.
        pixInd = sub2ind([numRows numCols],line0(:,1),line0(:,2));
        pixVals = voteMat(pixInd);
        pcount = numel(find(pixVals>vote(j)));   

        if(pcount>0) 
            voteMat(r(j),c(j)) = 0;         % NMS        
        end
 
    elseif(orientation==135)
        % get the list of all pixels in the bar
        % 1: the line through (c,r)    
        b = c+r;        % b is the y-intercept of the straight line y=ax+b going through (c,r)
        cStart = c + ceil(cosd(orientation)*barLength/2);
        cEnd = c - ceil(cosd(orientation)*barLength/2);

        cRange = cStart:cEnd;
        rRange = b - cRange;    % from the eqn of the straight line

        line0 = [rRange' cRange'];
        % for each additional line
        nLinePairs = (barWidth-1)/2; % number of additional pairs of lines
        for i = 1:nLinePairs
            rUp = r-i;
            bUp = c+rUp;
            cStartUp = cStart-i;
            cEndUp = cEnd-i;
            cRangeUp = cStartUp:cEndUp;
            rRangeUp = bUp-cRangeUp;

            rDown = r+i;
            bDown = c+rDown;
            cStartDown = cStart+i;
            cEndDown = cEnd+i;
            cRangeDown = cStartDown:cEndDown;
            rRangeDown = bDown-cRangeDown;

            lineUp = [rRangeUp' cRangeUp'];
            lineDown = [rRangeDown' cRangeDown'];

            line0 = [line0; lineUp; lineDown];        
        end
        % line0 contains (r,c) of each point in the bar as a column vector
        % count the pixels in line0 above the vote of this pixel
        pixInd = sub2ind([numRows numCols],line0(:,1),line0(:,2));
        pixVals = voteMat(pixInd);
        pcount = numel(find(pixVals>vote(j)));   

        if(pcount>0) 
            voteMat(r(j),c(j)) = 0;         % NMS        
        end
        
    end
        
end
