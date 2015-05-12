function dist = getRelativeDistance(curve,startVal,maxVal,value)

dist = interp1(curve,startVal:maxVal,value);