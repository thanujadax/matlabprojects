function h = plotLocalHoughSpaces(img,localHoughSpaces,T,R,maxLocalLines,...
                    thresholdFraction,houghSupNHood,fillGap,minLength,...
                    bb,slidingDist)

[rows cols] = size(localHoughSpaces);
idx = 1;
h = figure(200);

for i=1:rows
    startRow = (i-1)*(bb-slidingDist) + 1;
    stopRow = startRow + bb -1;
    for j=1:cols
        startCol = (j-1)*(bb-slidingDist) + 1;
        stopCol = startCol + bb -1;
        imgPatch = img(startRow:stopRow,startCol:stopCol);
        % extract lines from this houghSpace i,j
        H = cell2mat(localHoughSpaces(i,j));
        P  = houghpeaks(H,maxLocalLines,...
            'threshold',thresholdFraction*max(max(H)),...
            'NHoodSize',houghSupNHood);
       lines = houghlines(imgPatch,T,R,P,'FillGap',fillGap,'MinLength',minLength); 
       % corresponding patch in real image       
       subplot(rows,cols,idx); 
       imagesc(imgPatch);
       hold on
       for k = 1:length(lines)
           xy = [lines(k).point1; lines(k).point2];
           plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

           % Plot beginnings and ends of lines
           plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
           plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
        end
        hold off
       
       idx = idx + 1;
    end
end

