function [borderGridCellInds,borderFaceInds,allFacesOfBorderCells]...
    = getBoundaryCellInds(numCellsR,numCellsC,numSections)

% Output:
%   borderGridCellInds - 
%   borderFaceInds - indexed in the sequential order of all faces 1-6 given
%   in the order of the corresponding gridCellInds. e.g. Face1 of cell2 is
%   indexed 7.

numCellsPerSection = numCellsR * numCellsC;

% section = ones(numCellsR,numCellsC);
cellIDsInSection = 1:numCellsPerSection;
[rSection,cSection] = ind2sub([numCellsR numCellsC],cellIDsInSection');

% All cells of section_1
sectionID = 1;
z_1 = ones(numCellsPerSection,1) .* sectionID;
inds_sec_1 = sub2ind([numCellsR numCellsC numSections],rSection,cSection,z_1);

% All cells of section_n
sectionID = numSections;
z_1 = ones(numCellsPerSection,1) .* sectionID;
inds_sec_n = sub2ind([numCellsR numCellsC numSections],rSection,cSection,z_1);

% get border cell IDs for a single section
[borderCellIDs_internalSection, borderFaceIDs_internalSection] ...
    = getBorderCellsForSingleSection(numCellsR,numCellsC);
numBorderCells_internalSection = numel(borderCellIDs_internalSection);

totNumBorderCells = (numSections-2) * numBorderCells_internalSection ...
                        + 2*numCellsPerSection;
totNumBorderFaces = numCellsR*numCellsC*2 ...
                    + numCellsR*numSections*2 ...
                    + numCellsC*numSections*2;

% init
borderGridCellInds = zeros(totNumBorderCells,1);
borderFaceInds = zeros(totNumBorderFaces,1);

% section 1
borderGridCellInds(1:numCellsPerSection) = inds_sec_1;
% section1: face1xy
startInd_face = 1;
stopInd_face = numCellsPerSection;
borderFaceInds(startInd_face:stopInd_face) = (inds_sec_1-1)*6 +1;

% sections 2 to n-1
cellIndStop = numCellsPerSection;
faceIndStop = numCellsPerSection;
for i=2:numSections-1
    % cell ID of borders
    cellIndStart = cellIndStop +1;
    cellIndStop = cellIndStop + numBorderCells_internalSection;
    borderCellIDs_i =  borderCellIDs_internalSection ...
                + numCellsPerSection*(i-1);
    borderGridCellInds(cellIndStart:cellIndStop) = borderCellIDs_i;
    
    borderFaceIDs_i = borderFaceIDs_internalSection + numCellsPerSection ...
                    + (i-2)*(numCellsR+numCellsC)*2;
    faceIndStart = faceIndStop + 1;
    faceIndStop = faceIndStop + numel(borderFaceIDs_internalSection);
    borderFaceInds(faceIndStart:faceIndStop) = borderFaceIDs_i;
        
end
% section n (end)
cellIndStart = cellIndStop +1;
cellIndStop = cellIndStop + numCellsPerSection;
borderGridCellInds(cellIndStart:cellIndStop) = inds_sec_n;

borderCellFacesMat = zeros(totNumBorderCells,6);
% face1
borderCellFacesMat(:,1) = (borderGridCellInds-1)*6 +1;
% face2
borderCellFacesMat(:,2) = (borderGridCellInds-1)*6 +2;
% face3
borderCellFacesMat(:,3) = (borderGridCellInds-1)*6 +3;
% face4
borderCellFacesMat(:,4) = (borderGridCellInds-1)*6 +4;
% face5
borderCellFacesMat(:,5) = (borderGridCellInds-1)*6 +5;
% face6
borderCellFacesMat(:,6) = (borderGridCellInds-1)*6 +6;

numBorderFaces = numel(borderGridCellInds) * 6;
allFacesOfBorderCells = reshape(borderCellFacesMat,numBorderFaces,1);

%% Supplementary functions

function [borderCells_section, borderFaceIDs_section] ...
                = getBorderCellsForSingleSection(numCellsR,numCellsC)

% to get the indices of border cells and faces relative to one section
            
% row 1
colIDseq = 1:numCellsC;
rowIDseq  = ones(1,numCellsC);
row_1_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');
% border faces for row 1: xz_face5
row1_face5_xz =(row_1_cellIDs-1)*6 +5;

% row n
rowIDseq  = ones(1,numCellsC) .* numCellsR;
row_n_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');
% border faces for row n: xz_face6
rown_face6_xz =(row_n_cellIDs-1)*6 +6;

% col 1
rowIDseq = 1:numCellsR;
colIDseq = ones(1,numCellsR);
col_1_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');
% border faces for col 1: yz_face3
col1_face3_yz =(col_1_cellIDs-1)*6 +3;

% col n
colIDseq = ones(1,numCellsR) .* numCellsC;
col_n_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');
% border faces for col n: yz_face4
coln_face4_yz =(col_n_cellIDs-1)*6 +4;

% union
borderCells_section = [row_1_cellIDs; row_n_cellIDs; col_1_cellIDs; col_n_cellIDs];
borderCells_section = unique(borderCells_section);
borderFaceIDs_section = [row1_face5_xz; rown_face6_xz; col1_face3_yz; coln_face4_yz];
