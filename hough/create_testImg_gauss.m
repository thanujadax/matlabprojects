% create test image consisting of line segments
numRows = 100;
numCols = 620;

barWidth = 9;
barLength = 29;
img = zeros(numRows,numCols);

bar = zeros(barWidth,barLength);
sigX = 30;
sigY = 4;
centerR = (barWidth+1)/2;
centerC = (barLength+1)/2;

orientations = 0:10:170;
colshift = 0;
rshift = 0;
c0 = floor((barLength+1)/2);
r0 = floor((barWidth+1)/2);

r = 24;
c = 25;

for i = 1:numel(orientations)
    orientation = orientations(i);
    % rotate
    gaussBar = gauss2d(bar,[sigX sigY],[centerC centerR]);
    if(orientation>0)
        orientation = 180 - orientation;
        gaussBar = imrotate(gaussBar,orientation);
    end
    
    [rows cols] = find(gaussBar>0.1);
    barInd = find(gaussBar>0.1);
    numel(rows)
    rShift = r - r0;
    cShift = c - c0;

    % shift and place in img
    rows = rows + rShift;
    cols = cols + cShift;
    pixelInd = sub2ind([numRows numCols],rows,cols);
    img(pixelInd) = gaussBar(barInd);
   % figure(110);imshow(bar)
    
   c = c + 33;
end

figure(22)
imshow(img)
colormap('gray')

imwrite(img,'testImgLines3.png','png')
