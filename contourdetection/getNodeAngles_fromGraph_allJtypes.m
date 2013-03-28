function jAnglesAll_alpha = getNodeAngles_fromGraph_allJtypes(junctionTypeListInds,...
    nodeInds,jEdgesAll,edges2pixels,sizeR,sizeC)
% returns a cell array.
% jAnglesAll{i} - each row corresponds to the set of angles for each
% junction of type i (type 1 = J2)

% Inputs:
%   ...
%   jEdgesAll - cell array with jEdgesAll{i} corresponding the edgeSet of
%   junction type i. each row in the edge set corresponds to a junction
%   instance of that type

[nre,nce] = size(edges2pixels);  % first column is the edgeID
edgepixels = edges2pixels(:,2:nce);

[maxNodesPerJtype, numJtypes] = size(junctionTypeListInds);
maxNumEdgesPerNode = numJtypes + 1;     % list starts with J2 (= two edges)

% jAnglesAll = zeros(maxNodesPerJtype,maxNumEdgesPerNode,numJtypes);
jAnglesAll_alpha = cell(1,numJtypes);

for dim=1:numJtypes
    if(jEdgesAll{dim}==0)
        jAnglesAll_alpha{dim} = 0;
    else
        jEdges0 = jEdgesAll{dim};
        % jEdges = jEdges0(jEdges0>0);
        [r,c] = find(jEdges0>0);
        rmax = max(r);
        cmax = max(c);
        jEdges = zeros(rmax,cmax);
        jEdges(r,c) = jEdges0(r,c);       
        [numJ,degree] = size(jEdges);
        jAngles = zeros(numJ,degree);
        for i=1:numJ
            % for each node
            edges_i = jEdges(i,:);
            nodeListInd = junctionTypeListInds(i,dim);% get the index of the node in concern
            nodeInd = nodeInds(nodeListInd); 
            [rNode,cNode] = ind2sub([sizeR sizeC],nodeInd);
            for j=1:degree
                % for each edge of this node
                edgeID = edges_i(j);
                if(edgeID~=0)
                    edgeListInd = find(edges2pixels(:,1)==edgeID);  
                    if(isempty(edgeListInd))
                        continue;
                    end
                    edgePixelInds0 = edgepixels(edgeListInd,:);
                    %edgePixelInds = edgePixelInds(edgePixelInds>0);
                    [r1,c1] = find(edgePixelInds0>0);
                    rmax = max(r1);
                    cmax = max(c1);
                    edgePixelInds = zeros(rmax,cmax);
                    edgePixelInds(r1,c1) = edgePixelInds0(r1,c1);
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
                    jAngles(i,j) = medianAlpha;
                end
                    
            end
        end
        % assign the jAngles for this Jtype into jAnglesAll(:,:,Jtype)
        jAnglesAll_alpha{dim} = jAngles;
    end
end