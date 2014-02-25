function [A,b,senseArray] = getILP3DGridConstraints(cellStats)

% version 1.0
% 2014.02.25

% Constraints for the 3D_ILP_GRID

% Inputs:
%   cellStats - number of grid cells along each dimension 
%       [numR,numC,numZ] = [numY,numX,numZ] 

% Outputs:
%   A - sparse matrix containing coefficients for ILP constraints
%   b - RHS for the constraints
%   senseArray - 

% Column structure of A (structure of variable vector
%   Each grid cell gives rise to 7 state variables
%   1 - cell internal state
%   2 - state of face x-y (front)
%   3 - state of face x-y (back)
%   4 - state of face y-z (left)
%   5 - state of face y-z (right)
%   6 - state of face x-z (top)
%   7 - state of face x-z (bottom)
%   

%% Init
numCells = cellStats(1) * cellStats(2) * cellStats(3);
numConstraints = numCells * 9;
% numVariables = numCells * 7;

b = zeros(numConstraints,1);
senseArray(1:numConstraints) = '=';

numNonZeros = 2*6*numCells + 4*3*numCells; % withough accounting for boundaries
% boundary cells give rise to a lower number of constraints there by lower
% number of nonZeros for A.
ii = zeros(numNonZeros,1);
jj = zeros(numNonZeros,1);
ss = zeros(numNonZeros,1);
%% Compatibility of cell state and face states, with the states of the neighbors
% A + a4 = D + d3
% the constraints are symmetrical. => on avg, 3 constraints per cube
% apply constraints only for faces 1,3,5 explicitly for each cell
% for face 2,4,6 the constraints would be added implicitly due to
% neighbors

% i - y (rows)
% j - x (cols)
% k - z sections
constraintCount = 0;
nzElementInd = 0;   % index for ii,jj,ss elements
for k=1:numZ
    for i=1:numR
        for j=1:numC
            cellInd = sub2ind([numR numC numC],i,j,k);
            
            % if face1 has a neighbor (before in z direction)
            % face 1
            if(k>0)
                face = 1;
                face_ind = getFaceIndAbsGivenCell(cellInd,face);
                face_neighborCellInd = sub2ind([numR numC numC],i,j,(k-1));
                neighborFace = 2;
                neighbor_face_ind = getFaceIndAbsGivenCell...
                        (face_neighborCellInd,neighborFace);
                    
                % A row ID (constraint number)
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0;
                % senseArray(constraintCount) = '='; % default
                
                % A columns
                % cellID + cellID_face1_ind = negh_cellID + ...
                % neighCellID_face2_ind
                
                % cellID
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = cellID;
                ss(nzElementInd) = 1; % coefficient for A
                
                % cellID_face1_ind
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = face_ind;
                ss(nzElementInd) = 1; % coefficient for A
                
                % neighbor_cellID
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = face_neighborCellInd;
                ss(nzElementInd) = -1; % coefficient for A
                
                % neighbor_face2_ind
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = neighbor_face_ind;
                ss(nzElementInd) = -1; % coefficient for A
                
            end
            
            % if face 3 has a neighbor (to the left)
            % face 3
            if(j>0)
                face = 3;
                face_ind = getFaceIndAbsGivenCell(cellInd,face);
                face_neighborCellInd = sub2ind([numR numC numC],i,(j-1),k);
                
                neighborFace = 4;
                neighbor_face_ind = getFaceIndAbsGivenCell...
                        (face_neighborCellInd,neighborFace);
                    
                % A row ID (constraint number)
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0;
                % senseArray(constraintCount) = '='; % default
                
                
                % A columns
                % cellID + cellID_face1_ind = negh_cellID + ...
                % neighCellID_face2_ind
                
                % cellID
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = cellID;
                ss(nzElementInd) = 1; % coefficient for A
                
                % cellID_face1_ind
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = face_ind;
                ss(nzElementInd) = 1; % coefficient for A
                
                % neighbor_cellID
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = face_neighborCellInd;
                ss(nzElementInd) = -1; % coefficient for A
                
                % neighbor_face2_ind
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = neighbor_face_ind;
                ss(nzElementInd) = -1; % coefficient for A
                
            end
            
            % if face 5 has a neighbor (above)
            % face 5
            if(i>0)
                face = 5;
                face_ind = getFaceIndAbsGivenCell(cellInd,face);
                face_neighborCellInd = sub2ind([numR numC numC],(i-1),j,k);
                
                neighborFace = 6;
                neighbor_face_ind = getFaceIndAbsGivenCell...
                        (face_neighborCellInd,neighborFace);
                    
                % A row ID (constraint number)
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0;
                % senseArray(constraintCount) = '='; % default
                
                % A columns
                % cellID + cellID_face1_ind = negh_cellID + ...
                % neighCellID_face2_ind
                
                % cellID
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = cellID;
                ss(nzElementInd) = 1; % coefficient for A
                
                % cellID_face1_ind
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = face_ind;
                ss(nzElementInd) = 1; % coefficient for A
                
                % neighbor_cellID
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = face_neighborCellInd;
                ss(nzElementInd) = -1; % coefficient for A
                
                % neighbor_face2_ind
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = neighbor_face_ind;
                ss(nzElementInd) = -1; % coefficient for A
                
                
            end            
        end
    end
end

%% Compatibility of cell state and it's own face states
% A + ai <= 1; 6 per cube
for k=1:numZ
    for i=1:numR
        for j=1:numC
            cellInd = sub2ind([numR numC numC],i,j,k);
            
            % face 1
            constraintCount = constraintCount + 1;
            b(constraintCount) = 1.1;
            senseArray(constraintCount) = '<';
            face = 1;
            face_ind = getFaceIndAbsGivenCell(cellInd,face);

            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = cellInd;
            ss(nzElementInd) = 1; % coefficient for A
            
            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = face_ind;
            ss(nzElementInd) = 1; % coefficient for A
            
            % face 2
            constraintCount = constraintCount + 1;
            b(constraintCount) = 1.1;
            senseArray(constraintCount) = '<';
            face = 2;
            face_ind = getFaceIndAbsGivenCell(cellInd,face);

            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = cellInd;
            ss(nzElementInd) = 1; % coefficient for A
            
            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = face_ind;
            ss(nzElementInd) = 1; % coefficient for A
            
            % face 3
            constraintCount = constraintCount + 1;
            b(constraintCount) = 1.1;
            senseArray(constraintCount) = '<';
            face = 3;
            face_ind = getFaceIndAbsGivenCell(cellInd,face);

            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = cellInd;
            ss(nzElementInd) = 1; % coefficient for A
            
            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = face_ind;
            ss(nzElementInd) = 1; % coefficient for A
            
            % face 4
            constraintCount = constraintCount + 1;
            b(constraintCount) = 1.1;
            senseArray(constraintCount) = '<';
            face = 4;
            face_ind = getFaceIndAbsGivenCell(cellInd,face);

            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = cellInd;
            ss(nzElementInd) = 1; % coefficient for A
            
            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = face_ind;
            ss(nzElementInd) = 1; % coefficient for A
            
            % face 5
            constraintCount = constraintCount + 1;
            b(constraintCount) = 1.1;
            senseArray(constraintCount) = '<';
            face = 5;
            face_ind = getFaceIndAbsGivenCell(cellInd,face);

            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = cellInd;
            ss(nzElementInd) = 1; % coefficient for A
            
            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = face_ind;
            ss(nzElementInd) = 1; % coefficient for A
            
            % face 6
            constraintCount = constraintCount + 1;
            b(constraintCount) = 1.1;
            senseArray(constraintCount) = '<';
            face = 6;
            face_ind = getFaceIndAbsGivenCell(cellInd,face);

            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = cellInd;
            ss(nzElementInd) = 1; % coefficient for A
            
            nzElementInd = nzElementInd + 1;
            ii(nzElementInd) = constraintCount;
            jj(nzElementInd) = face_ind;
            ss(nzElementInd) = 1; % coefficient for A
            
            
        end
    end
end

%% Create sparse output matrix
% ii - list of row indices
% jj - list of column indices
% ss - list of values 
% A(ii(k),jj(k)) = ss(k);
A = sparse(ii,jj,ss);

%% Supplementary functions
function faceInd = getFaceIndAbsGivenCell(cellInd,face)
% first variable of each set of 7 vars per cube is the cube internal state
faceInd = (cellInd-1)*7 + (face+1);