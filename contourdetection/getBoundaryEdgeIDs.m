function boundaryEdges = getBoundaryEdgeIDs(ws,edges2pixels)

% input: ws

% output:
%   boundaryEdges - edgeIDs. 

% the WS id of the boundary == 1

boundaryPix = (ws==1);

edgePix = getEdgePixForWsFace(boundaryPix,ws);
edgeListInds = getEdgeIDsFromPixels(edgePix,edges2pixels);

edgeIDList = edges2pixels(:,1);

boundaryEdges = edgeIDList(edgeListInds);