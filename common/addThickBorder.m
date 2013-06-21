function outImg = addThickBorder(inImg,margin,val)
% adds a thick border of pixel value val of size 'margin' to the image
% input image is in matrix form

[sizeR,sizeC] = size(inImg);

newR = sizeR + margin * 2;
newC = sizeC + margin * 2;

outImg = ones(newR,newC).*val;

% set the input image inside the outImg
startRow = margin + 1;
endRow = newR - margin;
startCol = margin + 1;
endCol = newC - margin;

outImg(startRow:endRow,startCol,endCol) = inImg;

