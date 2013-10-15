function [faceAdj,edges2regions,setOfRegionsMat,twoRegionEdgeIDs,wsIDsForRegions]...
    = getFaceAdjFromWS(ws,edges2pixels,b_imWithBorder)
% For all the faces in WS, get the corresponding regionID.
% RegionID is the sequence number of each region in setOfRegions

% Inputs:
% ws - 
% edges2pixels,
% b_imWithBorder - boolean var. 1, if image has a thick border

% Outputs: 
%   faceAdj - adjacency graph of the faces of the input planar graph.
%   values are edgeIDs
%   edges2cells - each row corresponds to an edge. the two columns give you
%   the two cells that are connected by that edge. The order of edges
%   present is according to listOfEdgeIDs
%   setOfCellsMat - each row corresponds to a cell and contains the set of
%   edges bounding that cell as a row vector with zero padding.
%   listOfEdgeIDs - the edges considered in the faceAdjMatrix. Each of these edges
%   connects a pair of cells

numWsFaces = max(max(ws));
% initialize
wsIDsForRegions = 1:numWsFaces;
setOfRegionsMat = [];

if(b_imWithBorder)
    start = 2; % the first ws region is the border. ignore it.
else
    start = 1;
end

for i=start:numWsFaces
    clear intPix_i edgePix_i edgeSet_i
    % get internal pixels
    intPix_i = (ws==i);
    % get edge pixels
    edgePix_i = getEdgePixForWsFace(intPix_i);
    % get edges (edgeIDs)
    edgeSet_i = getEdgeSetFromEdgePixSet(edgePix_i,edges2pixels);
    setOfRegionsMat = [setOfRegionsMat; edgeSet_i]; 
    % get regionID
end


% add the cellList (index) as the first col of setOfCells. This is done so
% that we can reuse getAdjacencyMat() to creage faceAdj.
setOfRegionsMat_2 = [wsIDsForRegions' setOfRegionsMat];
[faceAdj,edges2regions,~,twoRegionEdgeIDs] = getAdjacencyMat(setOfRegionsMat_2);


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



