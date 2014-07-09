function edgePix = getEdgePixForWsFace(intPix,ws)
% inputs
% intPix - logical pix inds of internal pixels of the relevant ws region

intPix_inds = find(intPix>0);
[sizeR,sizeC] = size(ws);
neighbors8_pixind = get8Neighbors(intPix_inds,sizeR,sizeC);
allNeighbors8List_pixind = unique(neighbors8_pixind);

wsVals = ws(allNeighbors8List_pixind);

% which pix inds are zero (ws edges)
edgePix = allNeighbors8List_pixind(wsVals==0);



