function regionIDs = getRegionIDsForAllWsFaces(ws,setOfRegions,edges2pixels,imWithBorder)
% For all the faces in WS, get the corresponding regionID.
% RegionID is the sequence number of each region in setOfRegions

% Inputs:
% setOfRegions: contains the list of EdgeIDs for each region identified

numWsFaces = max(max(ws));

if(imWithBorder)
    start = 2;
else
    start = 1;
end

for i=start:numWsFaces
    % get internal pixels
    intPix_i = (ws==i);
    % get edge pixels
    edgePix_i = getEdgePixForWsFace(intPix_i);
    % get edges (edgeIDs)
    edgeSet_i = getEdgeSetFromEdgePixSet(edgePix_i,edges2pixels);
    
    % get regionID
end


function edgePix = getEdgePixForWsFace(intPix)

[sizeR, sizeC] = size(ws);
[r,c] = find(intPix);
topRow = min(r);
botRow = max(r);
leftCol = min(c);
rightCol = max(c);

% each row
for j=topRow:botRow
    row_j = find(intPix(j,:));
    lc_j = min(row_j);
    rc_j = max(row_j);
    el = lc_j -1;
    if(el>0)
        % it's an edge
        edgePix = [edgePix; row_j el];
    end
    er = rc_j +1;
    if(er<=sizeC)
        % it's an edge
        edgePix = [edgePix; row_j er];
    end
end

% each col
for j=leftCol:rightCol
    col_j = find(indPix(:,j));
    tr_j = min(col_j);
    br_j = max(col_j);
    et = tr_j -1;
    if(et>0)
        % it's an edge
        edgePix = [edgePix; et col_j];
    end
    eb = br_j +1;
    if(eb<=sizeR)
        % it's an edge
        edgePix = [edgePix; eb col_j];
    end
end

edgePix = sub2ind([sizeR sizeC],edgePix(:,1),edgePix(:,2));
edgePix = unique(edgePix);

function edgeSet = getEdgeSetFromEdgePixSet(edgePix,edges2pixels)
numEdgePix = numel(edgePix);
edgepixels = edges2pixels;
edgepixels(:,1) = [];
[x1,x2] = ismember(edgepixels,edgePix);
[~,x1sum] = sum(x1,2); % sum each row
edgeSet = edges2pixels((x1sum>0),1);



