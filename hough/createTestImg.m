% create test image consisting of line segments
rows = 48;
cols = 48;

img = zeros(rows,cols);

startRow = 1;
stopRow = 24;
startCol = 1;
stopCol = 24;

% block 1: vertical line in the center
img(25:48,12:13) = 1;

% block 2: horizontal line in the center
img(36:37,25:48) = 1;

% block 3: diagonal from bottom left
j=1;
for i=1:24
    img(i,j)=1;
    if(j>1)
        img(i,j-1)=1;
    end
    if(j<24)
        img(i,j+1)=1;
    end
    j=j+1;
end

% block 4: diagonal from bottom right
j=25;
for i=24:-1:1
    img(i,j)=1;
    if(j>1)
        img(i,j-1)=1;
    end
    if(j<48)
        img(i,j+1)=1;
    end
    j=j+1;
end

figure(22)
imagesc(img)
colormap('gray')

imwrite(img,'testImgLines.png','png')

