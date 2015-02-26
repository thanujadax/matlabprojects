function jAnglesAll_alpha = getNodeAngles_fromGraph_allJtypes(junctionTypeListInds,...
    nodeInds,jEdgesAll,edges2pixels,sizeR,sizeC,edges2nodes,connectedJunctionIDs)
% returns a cell array.
% jAnglesAll{i} - each row corresponds to the set of angles for each
% junction of type i (type 1 = J2)

% Inputs:
%   ...
%   jEdgesAll - cell array with jEdgesAll{i} corresponding the edgeSet of
%   junction type i. each row in the edge set corresponds to a junction
%   instance of that type
MAX_NUM_PIXELS = 10;  % maximum number of pixels away from the node used to calculate alpha
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
        jEdges = jEdgesAll{dim};
        % jEdges = jEdges0(jEdges0>0);
        % [r,c] = find(jEdges0>0);
        % rmax = max(r);
        % cmax = max(c);
        % jEdges = zeros(rmax,cmax);
        % jEdges(r,c) = jEdges0(r,c);       
        [numJ,degree] = size(jEdges);
        jAngles = zeros(numJ,degree);
        for i=1:numJ
            % for each node
            edges_i = jEdges(i,:);
            % start debug
%             if(i==42)
%                 aa = 99;
%             end
            % end debug
            nodeListInd = junctionTypeListInds(i,dim);% get the index of the node in concern
            nodeInd = nodeInds(nodeListInd); 
            
            [rNode,cNode] = ind2sub([sizeR sizeC],nodeInd);
            for j=1:degree
            % TODO: If the node is a cluster node determine the angle based on
            % the neares pixel of the clusternode to the edge in concern
            
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
                    rmax = max(r1); % number of pixels
                    cmax = max(c1);
                    edgePixelInds = zeros(rmax,cmax);
                    edgePixelInds(r1,c1) = edgePixelInds0(r1,c1);
                    % get the edge pixels(3) which are closest to the node i
                    nodePixels = getNodeEdgePixel(nodeInd,edgePixelInds,sizeR,sizeC,...
                                    MAX_NUM_PIXELS);
                    % get their orientation
                    [rP,cP] = ind2sub([sizeR sizeC],nodePixels');
                    numEdgePix = numel(nodePixels);
%                     orientations = zeros(numEdgePix,1);
                    if(numEdgePix<8)
                        % if the edge is not very long, only look at the pixels
                        % close to this node to determine the direction
                        % get the 2 junction nodes
                        edgeNodes = edges2nodes(edgeListInd,:);
                        if(edgeNodes(1)==nodeListInd)
                            node2ListInd = edgeNodes(2);
                        elseif(edgeNodes(2)==nodeListInd)
                            node2ListInd = edgeNodes(1);
                        else
                            disp('ERROR: getNodeAngles_fromGraph_allJtypes. node mismatch');
                        end
                        % calculate alpha based on these 2
                        nodeInd2 = nodeInds(node2ListInd); 
                        [rNode2,cNode2] = ind2sub([sizeR sizeC],nodeInd2);
                        if(isClusterNode(nodeInd,connectedJunctionIDs))
                            % is cluster node. pick the closest cluster
                            % pixel
                            % get all cluster pixels
                            clusterPixInds = getClusterPixInds(nodeInd,connectedJunctionIDs);
                            closestPixelInd = getClosestPixel(nodeInd2,clusterPixInds,sizeR,sizeC);
                            [rNode,cNode] = ind2sub([sizeR sizeC],closestPixelInd);
                        end
                        
                        y = rNode2 - rNode;
                        x = cNode2 - cNode;
                    else
                        % get alpha wrt the end edge pixels
                        % edge is too long. Take the direction wrt pixel at
                        % MAX_NUM_PIXELS away from current node which is set above
                        rp1 = rP(1);
                        rp2 = rP(end);
                        cp1 = cP(1);
                        cp2 = cP(end);
                        y = rp2 - rp1;
                        x = cp2 - cp1;
                    end
                    alpha = atan2d(y,x);
                    if(alpha<0)
                        alpha = alpha + 360;
                    end

                    jAngles(i,j) = alpha;
                end
                    
            end % for this node - all edges
            % check if there are duplicates. if yes, recalculate all the
            % angles fo that node wrt the node instead of the closest edge
            % pixel
            alphas = unique(jAngles(i,:));
            if(numel(alphas)<numel(jAngles(i,:)))
                % recalculate alphas
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
                        nodePixels = getNodeEdgePixel(nodeInd,edgePixelInds,sizeR,sizeC,...
                                        MAX_NUM_PIXELS);
                        % get their orientation
                        [rP,cP] = ind2sub([sizeR sizeC],nodePixels');
                        numEdgePix = numel(nodePixels);
    %                     orientations = zeros(numEdgePix,1);
                        if(numEdgePix==1)
                            % just one edge pixel
                            % get the 2 junction nodes
                            edgeNodes = edges2nodes(edgeListInd,:);
                            if(edgeNodes(1)==nodeListInd)
                                node2ListInd = edgeNodes(2);
                            elseif(edgeNodes(2)==nodeListInd)
                                node2ListInd = edgeNodes(1);
                            else
                                disp('ERROR: getNodeAngles_fromGraph_allJtypes. node mismatch');
                            end
                            % calculate alpha based on these 2
                            nodeInd2 = nodeInds(node2ListInd); 
                            [rNode2,cNode2] = ind2sub([sizeR sizeC],nodeInd2);
                            y = rNode2 - rNode;
                            x = cNode2 - cNode;
                        else
                            % get alpha wrt the node pixel
                            rp1 = rNode;
                            rp2 = rP(end);
                            cp1 = cNode;
                            cp2 = cP(end);
                            y = rp2 - rp1;
                            x = cp2 - cp1;
                        end
                        alpha = atan2d(y,x);
                        if(alpha<0)
                            alpha = alpha + 360;
                        end

                        jAngles(i,j) = alpha;
                    end

                end
                
                
                
            end
            
        end
        % assign the jAngles for this Jtype into jAnglesAll(:,:,Jtype)
        jAnglesAll_alpha{dim} = jAngles;
    end
end