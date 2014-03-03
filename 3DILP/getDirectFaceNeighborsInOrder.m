function setOfDirectNeighbors_6 = getDirectFaceNeighborsInOrder...
                                (thisCellID,thisCellSectionID,cellR,cellC,...
                                numCellsR,numCellsC,numSections)
                                
% Inputs:
%   thisCellID
%   thisCellSectionID
%   cellR - rowID of this cell
%   cellC - colID of this cell

% computes the cellLIDs of the direct neighbors for thisCellID in the order
% xy1(front),xy2(back),yz1(left),yz2(right),xz1(top),xz2(bottom)

setOfDirectNeighbors_6 = zeros(1,6);

% xy1(front)
face_z = thisCellSectionID -1;
face_r = cellR;
face_c = cellC;
setOfDirectNeighbors_6(1) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);

% xy2(back)
face_z = thisCellSectionID +1;
face_r = cellR;
face_c = cellC;
setOfDirectNeighbors_6(2) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);

% yz1(left)
face_z = thisCellSectionID;
face_r = cellR;
face_c = cellC -1;
setOfDirectNeighbors_6(3) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);

% yz2(right)
face_z = thisCellSectionID;
face_r = cellR;
face_c = cellC +1;
setOfDirectNeighbors_6(4) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);

% xz1(top)
face_z = thisCellSectionID;
face_r = cellR -1;
face_c = cellC;
setOfDirectNeighbors_6(5) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);

% xz2(bottom)
face_z = thisCellSectionID;
face_r = cellR +1;
face_c = cellC;
setOfDirectNeighbors_6(6) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
