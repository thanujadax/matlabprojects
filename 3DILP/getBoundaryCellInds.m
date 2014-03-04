function borderGridCellInds = getBoundaryCellInds(numCellsR,numCellsC,numSections)

% Output:
%   boundaryGridCellInds - 

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
borderCellIDs_internalSection = getBorderCellsForSingleSection(numCellsR,numCellsC);
numBorderCells_internalSection = numel(borderCellIDs_internalSection);

totNumBorderCells = (numSections-2) * numBorderCells_internalSection ...
                        + 2*numCellsPerSection;

% init
borderGridCellInds = zeros(totNumBorderCells,1);
% section 1
borderGridCellInds(1:numCellsPerSection) = inds_sec_1;
% sections 2 to n-1
cellIndStop = numCellsPerSection;
for i=2:numSections-1
    % cell ID of borders
    cellIndStart = cellIndStop +1;
    cellIndStop = cellIndStop + numBorderCells_internalSection;
    borderCellIDs =  borderCellIDs_internalSection ...
                + numCellsPerSection + numBorderCells_internalSection*(i-2);
    borderGridCellInds(cellIndStart:cellIndStop) = borderCellIDs;
end
% section n (end)
cellIndStart = cellIndStop +1;
cellIndStop = cellIndStop + numCellsPerSection;
borderGridCellInds(cellIndStart:cellIndStop) = inds_sec_n;

%% Supplementary functions

function borderCells_section = getBorderCellsForSingleSection(numCellsR,numCellsC)
% row 1
colIDseq = 1:numCellsC;
rowIDseq  = ones(1,numCellsC);
row_1_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');

% row n
rowIDseq  = ones(1,numCellsC) .* numCellsR;
row_n_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');

% col 1
rowIDseq = 1:numCellsR;
colIDseq = ones(1,numCellsR);
col_1_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');

% col n
colIDseq = ones(1,numCellsR) .* numCellsC;
col_n_cellIDs = sub2ind([numCellsR numCellsC],rowIDseq',colIDseq');

% union
borderCells_section = [row_1_cellIDs; row_n_cellIDs; col_1_cellIDs; col_n_cellIDs];
borderCells_section = unique(borderCells_section);
