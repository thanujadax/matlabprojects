function wsOut = assignRandomIndicesToWatershedTransform(wsIn)

maxInd = max(max(wsIn));
randomIndicesList = randperm(maxInd);

wsOut = wsIn;

% ws 1 will remain the same
for i=2:maxInd
    wsOut(wsIn==i) = randomIndicesList(i);
end