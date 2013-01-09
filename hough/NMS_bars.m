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

[numRows numCols] = size(voteMat);

totEl = numRows*numCols;
listEl = 1:totEl;
matIndexed = reshape(listEl,numRows,numCols);

for j=1:numPixels
    if(voteMat(r(j),c(j))==0)
        continue;
    end
    % get neighbor hood of (r(i),c(i))
    if(orientation==0)
        % horizontal bar
        if(voteMat(r(j),c(j))==0)
            continue;
        end
        startRow = r(j) - ceil(barWidth/2);
        stopRow = startRow + barWidth;
        startCol = c(j) - ceil(barLength/2);
        stopCol = startCol + barLength;

        pxls = find(voteMat>vote(j)); % pixels having a higher vote    

 
        allowedInd = matIndexed(startRow:stopRow,startCol:stopCol);
        numAllowed = numel(allowedInd);
        allowedIndList = reshape(allowedInd,numAllowed,1);
        
        pixInd = intersect(pxls,allowedIndList);
        
        if(numel(pixInd)>0)
            % this point is not a maximum. set it to zero
            voteMat(r(j),c(j)) = 0;
        else
            % if this pixel is a maximum and,
            % if there're any pixels having an equal vote inside the
            % neighborhood, make them zero    
            pxls = find(voteMat==vote(j));
            %[row col] = ind2sub([numRows numCols],pxls);            
            pixInd = intersect(pxls,allowedIndList);
            % remove the current pixel from this list
            thisInd = sub2ind([numRows numCols],r(j),c(j));
            pixInd = pixInd(pixInd~=thisInd);
            % pixInd = intersect(rowInd,colInd);
            if(numel(pixInd)>0)
                voteMat(pixInd)=0;
            end
        end
        
    elseif(orientation==90)
        % vertical bar

        startRow = r(j) - ceil(barLength/2);
        stopRow = startRow + barLength;
        startCol = c(j) - ceil(barWidth/2);
        stopCol = startCol + barWidth;

        pxls = find(voteMat>vote(j)); % pixels having a higher vote    

 
        allowedInd = matIndexed(startRow:stopRow,startCol:stopCol);
        numAllowed = numel(allowedInd);
        allowedIndList = reshape(allowedInd,numAllowed,1);
        
        pixInd = intersect(pxls,allowedIndList);
        
        if(numel(pixInd)>0)
            % this point is not a maximum. set it to zero
            voteMat(r(j),c(j)) = 0;
        else
            % if this pixel is a maximum and,
            % if there're any pixels having an equal vote inside the
            % neighborhood, make them zero    
            pxls = find(voteMat==vote(j));
            %[row col] = ind2sub([numRows numCols],pxls);            
            pixInd = intersect(pxls,allowedIndList);
            % remove the current pixel from this list
            thisInd = sub2ind([numRows numCols],r(j),c(j));
            pixInd = pixInd(pixInd~=thisInd);
            % pixInd = intersect(rowInd,colInd);
            if(numel(pixInd)>0)
                voteMat(pixInd)=0;
            end
        end   
        
    elseif(orientation==45)
        % get the list of all pixels in the bar
        % 1: the line through (c,r)    
        b = r(j)-c(j);        % b is the y-intercept of the straight line y=ax+b going through (c,r)
        cStart = c(j) - floor(cosd(orientation)*barLength/2);
        cEnd = c(j) + floor(cosd(orientation)*barLength/2);

        cRange = cStart:cEnd;
        rRange = cRange + b;    % from the eqn of the straight line

        line0 = [rRange' cRange'];
        % for each additional line
        nLinePairs = (barWidth-1)/2; % number of additional pairs of lines
        for i = 1:nLinePairs
            rUp = r(j)-i;
            bUp = rUp-c(j);
            cStartUp = cStart+i;
            cEndUp = cEnd+i;
            cRangeUp = cStartUp:cEndUp;
            rRangeUp = cRangeUp + bUp;

            rDown = r(j)+i;
            bDown = rDown-c(j);
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
        else
            % is a maximum
            % if this point is a maximum but there are other pixels in the
            % neighborhood with non-zero votes, make them zero.
            thisInd = sub2ind([numRows numCols],r(j),c(j));
            pixInd = pixInd(pixInd~=thisInd);
            voteMat(pixInd)=0;     % suppress all other neighbors

        end
 
    elseif(orientation==135)
        % get the list of all pixels in the bar
        % 1: the line through (c,r)    
        b = c(j)+r(j);        % b is the y-intercept of the straight line y=ax+b going through (c,r)
        cStart = c(j) + ceil(cosd(orientation)*barLength/2);
        cEnd = c(j) - ceil(cosd(orientation)*barLength/2);

        cRange = cStart:cEnd;
        rRange = b - cRange;    % from the eqn of the straight line

        line0 = [rRange' cRange'];
        % for each additional line
        nLinePairs = (barWidth-1)/2; % number of additional pairs of lines
        for i = 1:nLinePairs
            rUp = r(j)-i;
            bUp = c(j)+rUp;
            cStartUp = cStart-i;
            cEndUp = cEnd-i;
            cRangeUp = cStartUp:cEndUp;
            rRangeUp = bUp-cRangeUp;

            rDown = r(j)+i;
            bDown = c(j)+rDown;
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
