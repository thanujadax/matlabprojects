function dist = getRelativeDistance(image1,image2,curve,maxVal,value)

dist = interp1(curve,1:maxVal,value);