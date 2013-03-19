function pixelInds = getPixSetFromEdgeIDset(edgeInds,edges2pixels)

pixelInds = edges2pixels(edgeInds,:);

pixelInds = pixelInds(pixelInds>0);
