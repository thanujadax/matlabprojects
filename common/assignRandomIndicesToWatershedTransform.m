function wsOut = assignRandomIndicesToWatershedTransform(wsIn)

maxInd = max(max(wsIn));
randomIndicesList = randperm(maxInd);

wsOut = wsIn;

for i=1:maxInd
    wsOut(wsIn==i) = randomIndicesList(i);
end