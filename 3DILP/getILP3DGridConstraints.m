function [A,b,senseArray] = getILP3DGridConstraints(cellStats,borderCellInds)

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

%% Constraint enable/disable
con1_AdjCellsAndFaces = 1;
con2_CellAndFace = 1;
con3_continuity = 1;
con4_borderCells = 0;
con5_inOut = 0; % TODO: certain structures are inside others

% debug constraints
seq = 1:80;
con6_enforce_0 = 1; % set of cellInds to be labeled cellInterior (0)
    inactiveCellInds_con6 = [22 27 38 43 54 58]; 
con7_enforce_1 = 1; % set of cellInds to be labeled cellExterior (1)
    activeCellInds_con7 = setdiff(seq,inactiveCellInds_con6); 

%% Init
numR = cellStats(1);
numC = cellStats(2);
numZ = cellStats(3);

numCells = cellStats(1) * cellStats(2) * cellStats(3);
numCellsPerSection = numR*numC;
numBorderCells = numel(borderCellInds);

% init: precalculate number of constraints
numRowsCon1 = 0;
numRowsCon2 = 0;
numRowsCon3 = 0;
numRowsCon4 = 0;
numRowsCon6 = 0;
numRowsCon7 = 0;

% init: number of nonzero elements for each contraint type
numNZCon1 = 0;
numNZCon2 = 0;
numNZCon3 = 0;
numNZCon4 = 0;
numNZCon6 = 0;
numNZCon7 = 0;

if(con1_AdjCellsAndFaces)
    numRowsCon1 = numCells*3 - numR*numC - numR*numZ - numC*numZ;    
    numNZperRow = 4;
    numNZCon1 = numRowsCon1 * numNZperRow;
    disp('Constraint1: cell and face of direct neighbors - ACTIVE')
else
    disp('Constraint1: cell and face of direct neighbors - INACTIVE')
end
if(con2_CellAndFace)
    numRowsCon2 = numCells * 3;
    numNZperRow = 2;
    numNZCon2 = numRowsCon2 * numNZperRow;
end
if(con3_continuity)
    numRowsCon3 = (numCells-numBorderCells)*6*2;
    numNZperRow = 4;
    numNZCon3 = numRowsCon3 * numNZperRow;
end
if(con4_borderCells)
    numRowsCon4 = 1;
    numNZperRow = numBorderCells +1;
    numNZCon4 = numRowsCon4 * numNZperRow;
end
if(con6_enforce_0)
   numRowsCon6 = numel(inactiveCellInds_con6); 
   numNZperRow = 1;
   numNZCon6 = numRowsCon6 * numNZperRow;
   disp('Constraint6: enforced cellInterior - ACTIVE')
end
if(con7_enforce_1)
   numRowsCon7 = numel(activeCellInds_con7); 
   numNZperRow = 1;
   numNZCon7 = numRowsCon7 * numNZperRow;
   disp('Constraint7: enforced cellExterior - ACTIVE')
end

numConstraints = numRowsCon1 + numRowsCon2 + numRowsCon3 + numRowsCon4 ...
                    + numRowsCon6 + numRowsCon7;
% numConstraints = numCells * 6 - 3*numBorderCells ;

numVariables = numCells * 7;

b = zeros(numConstraints,1);
senseArray(1:numConstraints) = '=';

% numNonZeros = 2*3*numCells + 4*3*numCells; % without accounting for boundaries
numNonZeros = numNZCon1 + numNZCon2 + numNZCon3 + numNZCon4 + numNZCon6 ...
                + numNZCon7; 

% boundary cells give rise to a lower number of constraints there by lower
% number of nonZeros for A.
ii = zeros(numNonZeros,1);
jj = zeros(numNonZeros,1);
ss = zeros(numNonZeros,1);

%% Compatibility of cell state and face states, with the states of the neighbors
% Constraint 1.
if(con1_AdjCellsAndFaces)
    % for Z =1, we have one constraint less (2 inst of 3) per cell. Similar for
    % X= 1 and Y = 1 planes.
    % number = totNumberOfCells * 3 - numCellsPerSection.
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
                cellInd = sub2ind([numR numC numZ],i,j,k);

                % if face1 has a neighbor (before in z direction)
                % face 1
                if(k>1)
                    face = 1;
                    face_ind = getFaceIndAbsGivenCell(cellInd,face);
                    face_neighborCellInd = sub2ind([numR numC numC],i,j,(k-1));
                    face_neighborVarInd = (face_neighborCellInd-1)*7 +1;% colID of A
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
                    ii(nzElementInd) = constraintCount; % rowID of A
                    jj(nzElementInd) = (cellInd-1)*7 +1;% colID of A
                    ss(nzElementInd) = 1;        % coefficient for A

                    % cellID_face1_ind
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = face_ind;
                    ss(nzElementInd) = 1; % coefficient for A

                    % neighbor_cellID
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = face_neighborVarInd;
                    ss(nzElementInd) = -1; % coefficient for A

                    % neighbor_face2_ind
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = neighbor_face_ind;
                    ss(nzElementInd) = -1; % coefficient for A

                end

                % if face 3 has a neighbor (to the left)
                % face 3
                if(j>1)
                    face = 3;
                    face_ind = getFaceIndAbsGivenCell(cellInd,face);
                    face_neighborCellInd = sub2ind([numR numC numZ],i,(j-1),k);
                    face_neighborVarInd = (face_neighborCellInd-1)*7 +1;% colID of A

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
                    jj(nzElementInd) = (cellInd-1)*7 +1;
                    ss(nzElementInd) = 1; % coefficient for A

                    % cellID_face1_ind
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = face_ind;
                    ss(nzElementInd) = 1; % coefficient for A

                    % neighbor_cellID
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = face_neighborVarInd;
                    ss(nzElementInd) = -1; % coefficient for A

                    % neighbor_face2_ind
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = neighbor_face_ind;
                    ss(nzElementInd) = -1; % coefficient for A

                end

                % if face 5 has a neighbor (above)
                % face 5
                if(i>1)
                    face = 5;
                    face_ind = getFaceIndAbsGivenCell(cellInd,face);
                    face_neighborCellInd = sub2ind([numR numC numZ],(i-1),j,k);
                    face_neighborVarInd = (face_neighborCellInd-1)*7 +1;% colID of A

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
                    jj(nzElementInd) = (cellInd-1)*7 +1;
                    ss(nzElementInd) = 1; % coefficient for A

                    % cellID_face1_ind
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = face_ind;
                    ss(nzElementInd) = 1; % coefficient for A

                    % neighbor_cellID
                    nzElementInd = nzElementInd + 1;
                    ii(nzElementInd) = constraintCount;
                    jj(nzElementInd) = face_neighborVarInd;
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
end

%% Compatibility of cell state and it's own face states
% TODO: CHECK: border cells to be handled differently? No.
% A + ai <= 1; 3 per cube
% do only for faces 1,3,5. 
% faces 2,4,6 will automatically be constrained due to the coupling
% introduced by the previous constraint.

if(con2_CellAndFace)
    % number: each cell has 3 such constraints.
    for k=1:numZ
        for i=1:numR
            for j=1:numC
                cellInd = sub2ind([numR numC numZ],i,j,k);

                % face 1
                constraintCount = constraintCount + 1;
                b(constraintCount) = 1.1;
                senseArray(constraintCount) = '<';
                face = 1;
                face_ind = getFaceIndAbsGivenCell(cellInd,face);

                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = (cellInd-1)*7 +1;
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
                jj(nzElementInd) = (cellInd-1)*7 +1;
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
                jj(nzElementInd) = (cellInd-1)*7 +1;
                ss(nzElementInd) = 1; % coefficient for A

                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = face_ind;
                ss(nzElementInd) = 1; % coefficient for A

            end
        end
    end
end
%% Continuity constraint
if(con3_continuity)
    % For each cell, for each face, there are 2 constraints.
    % 2 since to have a closed volume each (2D) face (=tile) should
    % continue in 2 directions (consistently across the entire system)
    
    % border cells are ignored
    
    for k=2:(numZ-1)
        for j=2:(numC-1)
            for i=2:(numR-1)
                    % thisCell = cell_y

                cellInd_y = sub2ind([numR numC numZ],i,j,k);
                cellInd_z = sub2ind([numR numC numZ],i,(j+1),k);
                cellInd_h = sub2ind([numR numC numZ],i,(j+1),(k-1));
                cellInd_s = sub2ind([numR numC numZ],(i+1),j,k);
                cellInd_g = sub2ind([numR numC numZ],i,j,(k-1));
                cellInd_l = sub2ind([numR numC numZ],(i+1),j,(k+1));
                cellInd_d = sub2ind([numR numC numZ],i,(j+1),(k+1));
                cellInd_c = sub2ind([numR numC numZ],i,j,(k+1));
                cellInd_k = sub2ind([numR numC numZ],i,(j-1),(k+1));
                cellInd_p = sub2ind([numR numC numZ],(i+1),(j-1),k);
                cellInd_t = sub2ind([numR numC numZ],(i+1),(j+1),k);
                cellInd_a = sub2ind([numR numC numZ],(i-1),j,(k+1));
                cellInd_x = sub2ind([numR numC numZ],(i-1),(j+1),k);
                cellInd_m = sub2ind([numR numC numZ],(i-1),(j-1),k);
                cellInd_v = sub2ind([numR numC numZ],(i+1),j,(k-1));
                cellInd_w = sub2ind([numR numC numZ],(i-1),j,k);
                cellInd_e = sub2ind([numR numC numZ],(i-1),j,(k-1));
                cellInd_n = sub2ind([numR numC numZ],i,(j-1),k);
                cellInd_r = sub2ind([numR numC numZ],i,(j-1),(k-1));

                faceInd_y1 = getFaceIndAbsGivenCell(cellInd_y,1);
                faceInd_y2 = getFaceIndAbsGivenCell(cellInd_y,2);
                faceInd_y3 = getFaceIndAbsGivenCell(cellInd_y,3);
                faceInd_y4 = getFaceIndAbsGivenCell(cellInd_y,4);
                faceInd_y5 = getFaceIndAbsGivenCell(cellInd_y,5);
                faceInd_y6 = getFaceIndAbsGivenCell(cellInd_y,6);

                faceInd_s2 = getFaceIndAbsGivenCell(cellInd_s,2);
                faceInd_s4 = getFaceIndAbsGivenCell(cellInd_s,4);

                faceInd_l5 = getFaceIndAbsGivenCell(cellInd_l,5);

                faceInd_z2 = getFaceIndAbsGivenCell(cellInd_z,2);
                faceInd_z5 = getFaceIndAbsGivenCell(cellInd_z,5);

                faceInd_d3 = getFaceIndAbsGivenCell(cellInd_d,3);

                faceInd_c3 = getFaceIndAbsGivenCell(cellInd_c,3);
                faceInd_c5 = getFaceIndAbsGivenCell(cellInd_c,5);

                faceInd_k1 = getFaceIndAbsGivenCell(cellInd_k,1);

                faceInd_p4 = getFaceIndAbsGivenCell(cellInd_p,4);

                faceInd_t5 = getFaceIndAbsGivenCell(cellInd_t,5);

                faceInd_a1 = getFaceIndAbsGivenCell(cellInd_a,1);

                faceInd_x3 = getFaceIndAbsGivenCell(cellInd_x,3);

                faceInd_h2 = getFaceIndAbsGivenCell(cellInd_h,2);

                faceInd_g4 = getFaceIndAbsGivenCell(cellInd_g,4);
                faceInd_g6 = getFaceIndAbsGivenCell(cellInd_g,6);

                faceInd_v2 = getFaceIndAbsGivenCell(cellInd_v,2);

                faceInd_m6 = getFaceIndAbsGivenCell(cellInd_m,6);

                faceInd_r4 = getFaceIndAbsGivenCell(cellInd_r,4);

                faceInd_n1 = getFaceIndAbsGivenCell(cellInd_n,1);
                faceInd_n6 = getFaceIndAbsGivenCell(cellInd_n,6);

                faceInd_e6 = getFaceIndAbsGivenCell(cellInd_e,6);

                faceInd_w1 = getFaceIndAbsGivenCell(cellInd_w,1);
                faceInd_w3 = getFaceIndAbsGivenCell(cellInd_w,3);
                                               
                %% face1 (this = y) constr1 -
                % this cell is y. 
                % y1 - w1 -y5 -e6 <= 0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';               
                % y1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y1;
                ss(nzElementInd) = 1; % coefficient for A
                % -w1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_w1;
                ss(nzElementInd) = -1; % coefficient for A
                % -y5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y5;
                ss(nzElementInd) = -1; % coefficient for A
                % -e6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_e6;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face1 (this = y) constr2 -
                % y1 - n1 - y3 - r4 <=0        
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';                
                % y1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y1;
                ss(nzElementInd) = 1; % coefficient for A
                % -n1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_n1;
                ss(nzElementInd) = -1; % coefficient for A
                % -y3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y3;
                ss(nzElementInd) = -1; % coefficient for A
                % -r4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_r4;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face2 (this = y) constr1
                % y2 -s2 -y6 -l5 <=0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y2;
                ss(nzElementInd) = 1; % coefficient for A
                % -s2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_s2;
                ss(nzElementInd) = -1; % coefficient for A
                % -y6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y6;
                ss(nzElementInd) = -1; % coefficient for A
                % -l5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_l5;
                ss(nzElementInd) = -1; % coefficient for A
                               
                %% face2 (this = y) constr2
                % y2 -z2 -y4 -d3 <=0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y2;
                ss(nzElementInd) = 1; % coefficient for A
                % -z2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_z2;
                ss(nzElementInd) = -1; % coefficient for A
                % -y4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y4;
                ss(nzElementInd) = -1; % coefficient for A
                % -d3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_d3;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face3 constr1
                % y3 -c3 - y2 - k1 <=0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y3;
                ss(nzElementInd) = 1; % coefficient for A
                % -c3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_c3;
                ss(nzElementInd) = -1; % coefficient for A
                % -y2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y2;
                ss(nzElementInd) = -1; % coefficient for A
                % -k1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_k1;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face3 constr2
                % y3 -w3 - y5 - m6 <=0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y3;
                ss(nzElementInd) = 1; % coefficient for A
                % -w3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_w3;
                ss(nzElementInd) = -1; % coefficient for A
                % -y5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y5;
                ss(nzElementInd) = -1; % coefficient for A
                % -m6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_m6;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face4 constr 1
                % y4 - g4 - y1 - h2 <= 0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y4;
                ss(nzElementInd) = 1; % coefficient for A
                % -g4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_g4;
                ss(nzElementInd) = -1; % coefficient for A
                % -y1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y1;
                ss(nzElementInd) = -1; % coefficient for A
                % -h2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_h2;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face4 constr 2
                % y4 - s4 - y6 - t5 <=0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y4;
                ss(nzElementInd) = 1; % coefficient for A
                % -s4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_s4;
                ss(nzElementInd) = -1; % coefficient for A
                % -y6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y6;
                ss(nzElementInd) = -1; % coefficient for A
                % -t5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_t5;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face5 constr1
                % y5 - c5 - y2 - a1 <= 0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y5;
                ss(nzElementInd) = 1; % coefficient for A
                % -c5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_c5;
                ss(nzElementInd) = -1; % coefficient for A
                % -y2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y2;
                ss(nzElementInd) = -1; % coefficient for A
                % -a1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_a1;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face5 constr2
                % y5 - z5 - y4 - x3 <= 0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y5;
                ss(nzElementInd) = 1; % coefficient for A
                % -z5
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_z5;
                ss(nzElementInd) = -1; % coefficient for A
                % -y4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y4;
                ss(nzElementInd) = -1; % coefficient for A
                % -x3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_x3;
                ss(nzElementInd) = -1; % coefficient for A
                
                %% face6 constr1
                % y6 - g6 - y1 - v2 <= 0
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y6;
                ss(nzElementInd) = 1; % coefficient for A
                % -g6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_g6;
                ss(nzElementInd) = -1; % coefficient for A
                % -y1
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y1;
                ss(nzElementInd) = -1; % coefficient for A
                % -v2
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_v2;
                ss(nzElementInd) = -1; % coefficient for A
                %% face6 constr2
                % y6 - n6 - y3 -p4
                constraintCount = constraintCount + 1;
                b(constraintCount) = 0.1;
                senseArray(constraintCount) = '<';
                % y6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y6;
                ss(nzElementInd) = 1; % coefficient for A
                % -n6
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_n6;
                ss(nzElementInd) = -1; % coefficient for A
                % -y3
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_y3;
                ss(nzElementInd) = -1; % coefficient for A
                % -p4
                nzElementInd = nzElementInd + 1;
                ii(nzElementInd) = constraintCount;
                jj(nzElementInd) = faceInd_p4;
                ss(nzElementInd) = -1; % coefficient for A
                
            end
        end
    end
    
end
%% border cell constraint
% the state of the boundary cells are fixed to '1'
if(con4_borderCells)
    % get the colIDs for border cells
end

%% enforce a set of cells to be inactive (=0, cellInterior)
if(con6_enforce_0)
   numInactiveCells = numel(inactiveCellInds_con6);
   for n=1:numInactiveCells       
        cellVariableInd = (inactiveCellInds_con6(n)-1)*7 +1;
        constraintCount = constraintCount + 1;
        b(constraintCount) = 0; % cell interior is 0
        % senseArray(constraintCount) = '='; % default '='
        nzElementInd = nzElementInd + 1;
        ii(nzElementInd) = constraintCount;
        jj(nzElementInd) = cellVariableInd;
        ss(nzElementInd) = 1;
   end
end

%% enforce a set of cells to be active (=1, cellExterior)
if(con7_enforce_1)
   numActiveCells = numel(activeCellInds_con7);
   for n=1:numActiveCells       
        cellVariableInd = (activeCellInds_con7(n)-1)*7 +1;
        constraintCount = constraintCount + 1;
        b(constraintCount) = 1; % cell exterior is 1
        % senseArray(constraintCount) = '='; % default '='
        nzElementInd = nzElementInd + 1;
        ii(nzElementInd) = constraintCount;
        jj(nzElementInd) = cellVariableInd;
        ss(nzElementInd) = 1;
   end
end
%% Create sparse output matrix
% ii - list of row indices
% jj - list of column indices
% ss - list of values 
% A(ii(k),jj(k)) = ss(k);
% remove the unwanted zeros at the end of each ii,jj, and ss.
ii = ii(ii>0);
jj = jj(jj>0);
numNonZeros = numel(ii);
ss = ss(1:numNonZeros);
A = sparse(ii,jj,ss,constraintCount,numVariables);
% if(numel(b)>constraintCount)
%     % trim b
%     b(constraintCount+1:end) = [];
%     % senseArray
% end
%% Supplementary functions
function faceInd = getFaceIndAbsGivenCell(cellInd,face)
% first variable of each set of 7 vars per cube is the cube internal state
faceInd = (cellInd-1)*7 + (face+1);