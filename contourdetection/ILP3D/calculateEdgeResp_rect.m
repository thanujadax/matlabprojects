function edgeResp = calculateEdgeResp_rect...
            (orientation,edgeID,OFR,eps_orientations,edges2pixels)
%   calculates edge response in a rectangular region centered around
%   the mid point of the edge for the given orientation +- epsilon,
%   (col1 of edgeOrientations)

% calculate the relevant pix inds, in the proximity of this edge
edgePixels_i = edges2pixels(edgeID,:);
edgePixels_i(1) = []; % first element is the edgeID
pixIndOfr = getEdgeCoverage_rect(edgePixels_i,sizeR,sizeC,theta);
% get the edge response from OFR
edgeResp = calculateWeightedEdgeResp(pixIndOfr,OFR,...
                orientation,eps_orientations);

