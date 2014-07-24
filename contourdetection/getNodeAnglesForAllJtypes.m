function jAnglesAll = getNodeAnglesForAllJtypes(junctionTypeListInds,...
    nodeInds,jEdgesAll,edges2pixels,orientedScoreSpace3D,sizeR,sizeC,angleStep)
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
jAnglesAll = cell(1,numJtypes);

for dim=1:numJtypes
    if(jEdgesAll{dim}==0)
        jAnglesAll{dim} = 0;
    else
        % jEdges0 = jEdgesAll{dim};
        % [r,c] = find(jEdges0>0);
        jEdges = jEdgesAll{dim};
        % rmax = max(r);
        % cmax = max(c);
        % jEdges = zeros(rmax,cmax);
        % jEdges(r,c) = jEdges0(r,c);       
        [numJ,degree] = size(jEdges);
        jAngles = zeros(numJ,degree);
        numOrientations = size(orientedScoreSpace3D,3);
        for i=1:numJ
            % for each node
%             if(i==42)
%                 abc = 11;
%             end
            edges_i = jEdges(i,:);
            nodeListInd = junctionTypeListInds(i,dim);% get the index of the node in concern
            nodeInd = nodeInds(nodeListInd); 
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
                    [r,c] = ind2sub([sizeR sizeC],nodePixels');
                    numEdgePix = numel(nodePixels);
                    orientations = zeros(numEdgePix,numOrientations);
                    orientationIndex = zeros(1,numEdgePix);
                    for k=1:numEdgePix
                        orientations(k,1:numOrientations) = orientedScoreSpace3D(r(k),c(k),:);
                        [~,orientationIndex(k)] = max(orientations(k,:));
                    end
            %         [~,orientationIndex] = max(orientedScoreSpace3D(r,c,:),3);
                    edgeAngles = (orientationIndex-1).*angleStep;
                    if(numel(edgeAngles)==2)
                        % get the first angle
                        avgEdgeAngle = edgeAngles(1);
                    else
                        % get median
                        avgEdgeAngle = median(edgeAngles);
                    end
        %             jAngles(i,j) = avgEdgeAngle;
                    jAngles(i,j) = avgEdgeAngle;
                end
                    
            end
        end
        % assign the jAngles for this Jtype into jAnglesAll(:,:,Jtype)
        jAnglesAll{dim} = jAngles;
    end
end