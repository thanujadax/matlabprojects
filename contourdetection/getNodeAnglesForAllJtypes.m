function jAnglesAll = getNodeAnglesForAllJtypes(junctionTypeListInds,...
    nodeInds,jEdgesAll,edges2pixels,orientedScoreSpace3D,sizeR,sizeC,angleStep)
% returns a 3D array.
% jAngles(:,:,1) - each row corresponds to the set of angles for each
% junction of type 1 (= J2)
[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);
maxNumEdgesPerNode = numJtypes + 1;     % list starts with J2 (= two edges)

jAnglesAll = zeros(maxNodesPerJtype,maxNumEdgesPerNode,numJtypes);

for dim=1:numJtypes
    jEdges0 = jEdgesAll(:,:,dim);
    jEdges = jEdges0(jEdges0>0);
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
%             jAngles(i,j) = avgEdgeAngle;
            jAnglesAll(i,j,dim) = avgEdgeAngle;
        end
    end
    % assign the jAngles for this Jtype into jAnglesAll(:,:,Jtype)
    
end