function outImg = addThickBorder(inImg,margin,val)
% adds a thick border of pixel value val of size 'margin' to the image
% input image is in matrix form

maxPixVal = max(max(inImg));
if(maxPixVal>1)
    coef = val * 255;
else
    coef = val;
end

outImg = padarray(inImg,[margin margin],coef);



