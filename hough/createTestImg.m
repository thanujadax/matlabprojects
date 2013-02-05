% create test image consisting of line segments
numRows = 420;
numCols = 700;
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
c = 12;

for i = 1:numel(orientations)
    orientation = orientations(i);
    % rotate
    if(orientation>0)
        orientation = 180 - orientation;
        bar = imrotate(bar,orientation);
    end

    [rows cols] = find(bar);
    rShift = r - r0;
    cShift = c - c0;

    % shift and place in img
    rows = rows + rShift
    cols = cols + cShift
    pixelInd = sub2ind([numRows numCols],rows,cols);
    img(pixelInd) = 1;
    
    c = c + 15;
end

figure(22)
imagesc(img)
colormap('gray')

% imwrite(img,'testImgLines2.png','png')

