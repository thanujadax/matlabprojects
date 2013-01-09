function barInd = getBar(numRows,numCols,peakInd,barLength,barWidth,orientation)

[r c] = ind2sub([numRows numCols],peakInd);     % row col of the peak for which the bar should be found
totEl = numRows*numCols;
listEl = 1:totEl;
matIndexed = reshape(listEl,numRows,numCols);

if(orientation==0)
    startRow = r - ceil(barWidth/2);
    stopRow = startRow + barWidth;
    startCol = c - ceil(barLength/2);
    stopCol = startCol + barLength;
    
    bar = matIndexed(startRow:stopRow,startCol:stopCol);
    barInd = reshape(bar,numel(bar),1);

elseif(orientation==90)
        startRow = r - ceil(barLength/2);
        stopRow = startRow + barLength;
        startCol = c - ceil(barWidth/2);
        stopCol = startCol + barWidth;    
        
        bar = matIndexed(startRow:stopRow,startCol:stopCol);
        barInd = reshape(bar,numel(bar),1);
    
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
    barInd = sub2ind([numRows numCols],line0(:,1),line0(:,2));
    
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
    barInd = sub2ind([numRows numCols],line0(:,1),line0(:,2));
end