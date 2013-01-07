function pcount = pixelCountInBar(img,r,c,orientation,barLength,barWidth,grayThresh)
pcount = 0; % init

if(orientation==0)
    % horizontal bar
    startRow = r - ceil(barWidth/2);
    stopRow = startRow + barWidth;
    startCol = c - ceil(barLength/2);
    stopCol = startCol + barLength;
    
    [pxlR, pxlC] = find(img>grayThresh); % non-zero pixels    
    rowIndInd = find(pxlR>startRow && pxlR<stopRow);    
    colIndInd = find(pxlC>startCol && pxlC<stopCol);    
    pixInd = intersect(rowIndInd,colIndInd);
    
    pcount = length(pixInd);    
    
elseif(orientation==90)
    % vertical bar
    startRow = r - ceil(barLength/2);
    stopRow = startRow + barLength;
    startCol = c - ceil(barWidth/2);
    stopCol = startCol + barWidth;
    
    [pxlR, pxlC] = find(img>grayThresh); % non-zero pixels    
    rowIndInd = find(pxlR>startRow && pxlR<stopRow);    
    colIndInd = find(pxlC>startCol && pxlC<stopCol);    
    pixInd = intersect(rowIndInd,colIndInd);
    
    pcount = length(pixInd);

elseif(orientation==45)
    rEnd = r + ceil(sind(45)*barLength/2);
    cEnd = c + ceil(cosd(45)*barLength/2);
    
    b = c-r;
    % get the list of all pixels in the bar
    
    
    
    
elseif(orientation==135)
    
end