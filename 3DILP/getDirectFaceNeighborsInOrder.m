function [setOfDirectNeighbors_6,setOfDirectNeighbors2D_4] ...
                = getDirectFaceNeighborsInOrder...
                                (thisCellID,numCellsR,numCellsC,numSections)
                                
% Inputs:
%   thisCellID
%   thisCellSectionID
%   cellR - rowID of this cell
%   cellC - colID of this cell

% computes the cellLIDs of the direct neighbors for thisCellID in the order
% xy1(front),xy2(back),yz1(left),yz2(right),xz1(top),xz2(bottom)

% if there is no direct neighbor for a particular face, then it returns
% zero for that particular face.

[cellR,cellC,thisCellSectionID] = ind2sub([numCellsR numCellsC numSections],thisCellID);

setOfDirectNeighbors_6 = zeros(1,6);
setOfDirectNeighbors2D_4 = zeros(1,4);

% xy1(front)
if(thisCellSectionID>1)
    face_z = thisCellSectionID -1;
    face_r = cellR;
    face_c = cellC;
    setOfDirectNeighbors_6(1) ...
            = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end

% xy2(back)
if(thisCellSectionID<numSections)
    face_z = thisCellSectionID +1;
    face_r = cellR;
    face_c = cellC;
    setOfDirectNeighbors_6(2) ...
            = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end

% yz1(left)
if(cellC>1)
    face_z = thisCellSectionID;
    face_r = cellR;
    face_c = cellC -1;
    neighID = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
    setOfDirectNeighbors_6(3) = neighID;
    setOfDirectNeighbors2D_4(1) = neighID;
end


% yz2(right)
if(cellC<numCellsC)
    face_z = thisCellSectionID;
    face_r = cellR;
    face_c = cellC +1;
    neighID = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
    setOfDirectNeighbors_6(4) = neighID;
    setOfDirectNeighbors2D_4(2) = neighID;
end

% xz1(top)
if(cellR>1)
    face_z = thisCellSectionID;
    face_r = cellR -1;
    face_c = cellC;
    neighID = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
    setOfDirectNeighbors_6(5) = neighID;
    setOfDirectNeighbors2D_4(3) = neighID;
end

% xz2(bottom)
if(cellR<numCellsR)
    face_z = thisCellSectionID;
    face_r = cellR +1;
    face_c = cellC;
    neighID = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
    setOfDirectNeighbors_6(6) = neighID;
    setOfDirectNeighbors2D_4(4) = neighID;
end