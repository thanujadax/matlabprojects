function closestPixelInd = getClosestPixel(refPixInd,pixIndList,sizeR,sizeC)

[refR,refC] = ind2sub([sizeR sizeC],refPixInd);

[pixListR, pixListC] = ind2sub([sizeR sizeC],pixIndList);

distance = (pixListR - refR).^2 + (pixListC - refC).^2;

[sortedDists, sortedInd] = sort(distance);

sortedPixInds = pixIndList(sortedInd);

closestPixelInd = sortedPixInds(1);