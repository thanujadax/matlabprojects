function connectedJunctionIDs = getClusteredJunctions(wsJ)
% Inputs:
%   wsJ - matrix containing 1 for junction pixels detected. others 0.

% Output:
%   connectedJunctionIDs - N-by-2 array. each cluster of junctions assigned
%   a unique label from 1 to n. n = number of clusters.


% extract edges with zero pixel length
% junction nodes (pixels) which are next to each other
% wsJ contains 1 for junction nodes, others zero.

[sizeR,sizeC] = size(wsJ);
ind4J = find(wsJ);

eightNH_J = zeros(sizeR,sizeC);
numJunctionPixels = numel(ind4J);  % watershed junction pixels (3 and 4 Js)
for i=1:numJunctionPixels
    % calculate n.o. 4 neighbors
    ind = ind4J(i);
    [r c] = ind2sub([sizeR sizeC],ind);
    nh = zeros(1,8);
    % N1
    if((r-1)>0)
        nh(1) = wsJ(r-1,c);
    end
    % N2
    if((r+1)<=sizeR)
        nh(2) = wsJ(r+1,c);
    end
    % N3
    if((c-1)>0)
        nh(3) = wsJ(r,c-1);
    end
    % N4
    if((c+1)<=sizeC)
        nh(4) = wsJ(r,c+1);
    end
    % N5
    if((r-1)>0 && (c-1)>0)
        nh(5) = wsJ(r-1,c-1);
    end
    % N6
    if((r-1)>0 && (c+1)<=sizeC)
        nh(6) = wsJ(r-1,c+1);
    end
    % N7
    if((r+1)<=sizeR && (c-1)>0)
        nh(7) = wsJ(r+1,c-1);
    end
    % N8
    if((r+1)<=sizeR && (c+1)<=sizeC)
        nh(8) = wsJ(r+1,c+1);
    end
    eightNH_J(ind) = sum(nh);
end

indJClusterPixels = find(eightNH_J>0);  % gets junction pixels which has neighbors
% make a N-by-2 array of such junction pixels having immediate neighboring
% junctions
numJunctionClustPix = numel(indJClusterPixels);
if(numJunctionClustPix==0)
    % no clustered junction nodes
    connectedJunctionIDs = 0;
else
    connectedJunctionIDs = zeros(numJunctionClustPix,2);
    connectedJunctionIDs(:,1) = indJClusterPixels;
    junctionLabel = 0;
    for i=1:numJunctionClustPix
        % look for the neighbors and give them the same label
        jLabelCurrent = connectedJunctionIDs(i,2);
        if(jLabelCurrent==0)
            junctionLabel = junctionLabel + 1;
            % assign label to this junction and to its neighbors and its
            % neighbors neighbors
            junctionInd = connectedJunctionIDs(i,1);
            connectedJunctionIDs = labelJunctionClusterNeighbors(junctionInd,connectedJunctionIDs,junctionLabel,...
                        sizeR,sizeC,eightNH_J);
        end
    end
    % connectedJunctionIDs contain the same ID for each node that is connected
    % together with zero length edges
    % nodeZeroEdges - store node - edge1,edge2 etc for these zero length edges
end