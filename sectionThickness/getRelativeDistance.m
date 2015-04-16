function dist = getRelativeDistance(curve,maxVal,value)

dist = interp1(curve,1:maxVal,value);