function connectedJunctionIDs = labelJunctionCluster4Neighbors(junctionInd,connectedJunctionIDs,junctionLabel,...
                sizeR,sizeC,eightNH_J)
% assign new label to this junction and all its neighbors and the
% neighbors' neighbors

i = find(connectedJunctionIDs(:,1)==junctionInd);
currentJlabel = connectedJunctionIDs(i,2);
if(currentJlabel==0)    
    connectedJunctionIDs(i,2) = junctionLabel;

    % now deal with its neighbors
    % ind = connectedJunctionIDs(i,1);
    [r,c] = ind2sub([sizeR sizeC],junctionInd);
    % N1
    if((r-1)>0)
        isJunction = eightNH_J((r-1),c);
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],(r-1),c);
%             listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
%             connectedJunctionIDs(listIndNeigh,2) = junctionID;
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
    % N2
    if((r+1)<=sizeR)
        isJunction = eightNH_J((r+1),c);
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],(r+1),c);
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
    % N3
    if((c-1)>0)
        isJunction = eightNH_J(r,(c-1));
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],r,(c-1));
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
    % N4
    if((c+1)<=sizeC)
        isJunction = eightNH_J(r,(c+1));
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],r,(c+1));
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
    % N5
    if((r-1)>0 && (c-1)>0)
        isJunction = eightNH_J((r-1),(c-1));
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],(r-1),(c-1));
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
    % N6
    if((r-1)>0 && (c+1)<=sizeC)
        isJunction = eightNH_J((r-1),(c+1));
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],(r-1),(c+1));
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
    % N7
    if((r+1)<=sizeR && (c-1)>0)
        isJunction = eightNH_J((r+1),(c-1));
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],(r+1),(c-1));
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
    % N8
    if((r+1)<=sizeR && (c+1)<=sizeC)
        isJunction = eightNH_J((r+1),(c+1));
        if(isJunction)
            indNeigh = sub2ind([sizeR sizeC],(r+1),(c+1));
            connectedJunctionIDs = labelJunctionCluster4Neighbors(indNeigh,connectedJunctionIDs,...
                junctionLabel,sizeR,sizeC,eightNH_J);
        end
    end
end