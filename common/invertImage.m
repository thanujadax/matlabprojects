function invertedImg = invertImage(imgIn)

% % scaling 0-255
% imgIn = int32(imgIn);
% imgIn = imgIn./(max(max(imgIn)));
% imgIn = imgIn .* 255;

if(max(max(imgIn))>1)
    maxVal = 255;
else
    maxVal = 1;
end

imgIn = double(imgIn);
invertedImg = imgIn - maxVal;
invertedImg = invertedImg .* (-1);
% invertedImg = uint8(invertedImg);
