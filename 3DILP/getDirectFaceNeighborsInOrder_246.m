function setOfDirectNeighbors_3 = getDirectFaceNeighborsInOrder_246...
                                (thisCellID,thisCellSectionID,cellR,cellC,...
                                numCellsR,numCellsC,numSections)
                                
% Inputs:
%   thisCellID
%   thisCellSectionID
%   cellR - rowID of this cell
%   cellC - colID of this cell

% computes the cellLIDs of the direct neighbors for thisCellID in the order
% 2-4-6 = xy2(back),yz2(right),xz2(bottom)

% if there is no direct neighbor for a particular face, then it returns
% zero for that particular face.

setOfDirectNeighbors_3 = zeros(1,3);

% xy2(back)
if(thisCellSectionID<numSections)
    face_z = thisCellSectionID +1;
    face_r = cellR;
    face_c = cellC;
    setOfDirectNeighbors_3(1) ...
            = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end

% yz2(right)
if(cellC<numCellsC)
    face_z = thisCellSectionID;
    face_r = cellR;
    face_c = cellC +1;
    setOfDirectNeighbors_3(1) ...
        = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end

% xz2(bottom)
if(cellR<numCellsR)
    face_z = thisCellSectionID;
    face_r = cellR +1;
    face_c = cellC;
    setOfDirectNeighbors_3(1) ...
            = sub2ind([numCellsR numCellsC numSections],face_r,face_c,face_z);
end