function outImg = removeBorder(inImg,dim)

% returns the output image removing the border of inImg so that outImg has
% dimensions dim
% the centers of inImg and outImg are the same

% inputs 
% inImg = original input image
% dim = [sizeR sizeC] % number of rows and cols

sizeOrig = size(inImg);

borderRows = floor((sizeOrig(1) - dim(1))/2);
borderCols = floor((sizeOrig(2) - dim(2))/2);

startRow = borderRows + 1;
endRow = sizeOrig(1) - borderRows;
startCol = borderCols + 1;
endCol = sizeOrig(2) - borderCols;

outImg = inImg(startRow:endRow,startCol:endCol);

