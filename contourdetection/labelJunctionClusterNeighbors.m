function connectedJunctionIDs = labelJunctionClusterNeighbors(i,connectedJunctionIDs,junctionID,...
                sizeR,sizeC,eightNH_J)
% assign new label to this junction and all its neighbors and the
% neighbors' neighbors

connectedJunctionIDs(i,2) = junctionID;

% now deal with its neighbors
ind = connectedJunctionIDs(i,1);
[r c] = ind2sub([sizeR sizeC],ind);
% N1
if((r-1)>0)
    isJunction = eightNH_J((r-1),c);
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],(r-1),c);
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end
% N2
if((r+1)<=sizeR)
    isJunction = eightNH_J((r+1),c);
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],(r+1),c);
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end
% N3
if((c-1)>0)
    isJunction = eightNH_J(r,(c-1));
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],r,(c-1));
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end
% N4
if((c+1)<=sizeC)
    isJunction = eightNH_J(r,(c+1));
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],r,(c+1));
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end
% N5
if((r-1)>0 && (c-1)>0)
    isJunction = eightNH_J((r-1),(c-1));
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],(r-1),(c-1));
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end
% N6
if((r-1)>0 && (c+1)<=sizeC)
    isJunction = eightNH_J((r-1),(c+1));
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],(r-1),(c+1));
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end
% N7
if((r+1)<=sizeR && (c-1)>0)
    isJunction = eightNH_J((r+1),(c-1));
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],(r+1),(c-1));
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end
% N8
if((r+1)<=sizeR && (c+1)<=sizeC)
    isJunction = eightNH_J((r+1),(c+1));
    if(isJunction)
        indNeigh = sub2ind([sizeR sizeC],(r+1),(c+1));
        listIndNeigh = find(connectedJunctionIDs(:,1)==indNeigh);
        connectedJunctionIDs(listIndNeigh,2) = junctionID;
    end
end