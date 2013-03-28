function jAngles_alpha = getNodeAngles_fromGraph(jInd,nodeInds,jEdges,edges2pixels,...
                    sizeR,sizeC)
% returns an N-by-n array of angles (alphas). N = number of nodes, n = number of
% edges per node. The order is determined by jInd which contains the
% indices of the nodes (junctions).
% Here, the angle of an edge at a node is the actual orientation of an edge
% based on its physical location in the graph obtained from the watershed segmentation. 

% Inputs:
%   jInd = list of node indices considered in order
%   jEdges = the edges connected to each node given in jInd
%   edges2pixels = gives the pixels contained in each edge. First col has
%   the edgeInd
%   orientedScoreSpace3D = oriented filter response for each pixel

% Output:
%   jAngles = angles corresponding to the edges of each node as row vectors
%       the angle is determined by the actual physical orientation of the
%       edge according to the graph obtained from the watershed
%       segmentation.

[numJ,degree] = size(jEdges);
jAngles_alpha = zeros(numJ,degree);
for i=1:numJ
    % for each node
    edges_i = jEdges(i,:);
    nodeInd = nodeInds(jInd(i));
    [rNode,cNode] = ind2sub([sizeR sizeC],nodeInd);
    for j=1:degree
        % for each edge of this node
        edgeID = edges_i(j);
        edgePixelInds = edges2pixels(edgeID,:);
        edgePixelInds = edgePixelInds(edgePixelInds>0);
        % get the edge pixels(3) which are closest to the node i
        nodePixels = getNodeEdgePixel(nodeInd,edgePixelInds,sizeR,sizeC);
        % get their orientation
        [rP,cP] = ind2sub([sizeR sizeC],nodePixels');
        numEdgePix = numel(nodePixels);
        orientations = zeros(numEdgePix,1);
        for k=1:numEdgePix
            y = rP(k) - rNode;
            x = cP(k) - cNode;
            alpha = atan2d(y,x);
            if(alpha<0)
                alpha = alpha + 360;
            end
            orientations(k) = alpha;
        end

        medianAlpha = median(orientations);
        jAngles_alpha(i,j) = medianAlpha;
    end
end