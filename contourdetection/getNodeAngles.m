function jAngles = getNodeAngles(jInd,nodeInds,jEdges,edges2pixels,orientedScoreSpace3D,...
                    sizeR,sizeC,angleStep)
% returns an N-by-n array of angles. N = number of nodes, n = number of
% edges per node. The order is determined by jInd which contains the
% indices of the nodes (junctions)

% Inputs:
%   jInd = list of node indices considered in order
%   jEdges = the edges connected to each node given in jInd
%   edges2pixels = gives the pixels contained in each edge. First col has
%   the edgeInd
%   orientedScoreSpace3D = oriented filter response for each pixel

% Output:
%   jAngles = angles corresponding to the edges of each node as row vectors
%       the angle is determined by the maximum response obtained by the
%       oriented filter bank

[numJ,degree] = size(jEdges);
jAngles = zeros(numJ,degree);
numOrientations = size(orientedScoreSpace3D,3);
for i=1:numJ
    % for each node
    edges_i = jEdges(i,:);
    nodeInd = nodeInds(jInd(i));
    for j=1:degree
        % for each edge of this node
        edgeID = edges_i(j);
        edgePixelInds = edges2pixels(edgeID,:);
        edgePixelInds = edgePixelInds(edgePixelInds>0);
        % get the edge pixels(3) which are closest to the node i
        nodePixels = getNodeEdgePixel(nodeInd,edgePixelInds,sizeR,sizeC);
        % get their orientation
        [r,c] = ind2sub([sizeR sizeC],nodePixels');
        numEdgePix = numel(nodePixels);
        orientations = zeros(numEdgePix,numOrientations);
        for k=1:numEdgePix
            orientations(k,1:numOrientations) = orientedScoreSpace3D(r(k),c(k),:);
            [~,orientationIndex(k)] = max(orientations(k,:));
        end
%         [~,orientationIndex] = max(orientedScoreSpace3D(r,c,:),3);
        edgeAngles = (orientationIndex-1).*angleStep;
        avgEdgeAngle = median(edgeAngles);
        jAngles(i,j) = avgEdgeAngle;
    end
end