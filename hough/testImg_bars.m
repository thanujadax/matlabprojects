% create test image consisting of line segments
numRows = 50;
numCols = 400;
barWidth = 3;
barLength = 15;
img = zeros(numRows,numCols);

bar = ones(barWidth,barLength);
orientations = 0:10:170;
colshift = 0;
rshift = 0;
c0 = floor((barLength+1)/2);
r0 = floor((barWidth+1)/2);

r = 24;
c = 25;

for i = 1:numel(orientations)
    orientation = orientations(i)
    % rotate
    bar = ones(barWidth,barLength);
    if(orientation>0)
        orientation = 180 - orientation;
        bar = imrotate(bar,orientation);
    end
    
    [rows cols] = find(bar>0);
    numel(rows)
    rShift = r - r0;
    cShift = c - c0;

    % shift and place in img
    rows = rows + rShift;
    cols = cols + cShift;
    pixelInd = sub2ind([numRows numCols],rows,cols);
    img(pixelInd) = 1;
   % figure(110);imshow(bar)
    
   c = c + 20;
end

figure(22)
imshow(img)
colormap('gray')

imwrite(img,'testImgLines3.png','png')
