function invertedImg = invertImage(imgIn)

% scaling 0-255
imgIn = int32(imgIn);
imgIn = imgIn./(max(max(imgIn)));
imgIn = imgIn .* 255;

invertedImg = imgIn - 255;
invertedImg = invertedImg .* (-1);
invertedImg = uint8(invertedImg);
