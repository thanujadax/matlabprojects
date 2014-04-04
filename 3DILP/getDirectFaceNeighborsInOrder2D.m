function setOfDirectNeighbors2D_4 = getDirectFaceNeighborsInOrder2D...
                    (thisCellID,numCellsR,numCellsC,numSections)
                
% Inputs:


% Computes the cellIDs of the direct neighbors for thisCellID on the 2D
% section, in the following order:
%   (1)yz1(left)
%   (2)yz2(right)
%   (3)xz1(top)
%   (4)xz2(bottom)


[cellR,cellC,thisCellSectionID] = ind2sub([numCellsR numCellsC numSections],thisCellID);

setOfDirectNeighbors2D_4 = zeros(1,4);

% yz1(left) (face3)
if(cellC>1)
    face_z = thisCellSectionID;
    face_r = cellR;
    face_c = cellC -1;
    setOfDirectNeighbors2D_4(1) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end

% yz2(right) (face4)
if(cellC<numCellsC)
    face_z = thisCellSectionID;
    face_r = cellR;
    face_c = cellC +1;
    setOfDirectNeighbors2D_4(2) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end

% xz1(top) (face5)
if(cellR>1)
    face_z = thisCellSectionID;
    face_r = cellR -1;
    face_c = cellC;
    setOfDirectNeighbors2D_4(3) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end

% xz2(bottom) (face6)
if(cellR<numCellsR)
    face_z = thisCellSectionID;
    face_r = cellR +1;
    face_c = cellC;
    setOfDirectNeighbors2D_4(4) ...
            = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end