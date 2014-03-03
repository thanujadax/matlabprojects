function computeFmGridFaces(pathToFm,boundaryGridCellInds,...
                gridCIDs_sectionIDs_rootPixIDsRel,numZ,numCellsY,numCellsX)

% computes features for each face of eah grid cell
% takes into account the features of the direct neighbor of the face as
% well as the other neighbors around it.

% FEATURES ARE NOT CALCULATED FOR GRID BOUNDARY FACES!

% version 1.0

% TODO: version 1.0 takes into account only the direct neighbor of the
% face. Extend it to consider the features of the neighbors as well.

% each grid cell has 6 faces
% should not consider the boundary cells

% Face order
%   xy1(front),xy2(back),yz1(left),yz2(right),xz1(top),xz2(bottom)

% for each cell
% get the neighbors in order (starting with the direct neighbor)
% calculate differential features (direct neighbor)
% combine all features to one vector 
% save all faces in one file

% __________________
% |       |         |
% |       |         |
% |  A   a|b    B   |
% |       |         |
% |_______|_________|

% feature vector for face a = {f(A-B),f(A),f(B)}

% path to save section features
saveFilePath = fullfile(pathToFm,subDir_cellFaceFm);
pathToGridCellInteriorFeatures = fullfile(pathToFm,subDir_Fm);

% check if subdirectory exists. Create if not
checkAndCreateSubDir(pathToFm,subDir_cellInteriorFm);

% featureFile.mat name structure: fm_slice_%d.mat

numGridCellsTot = size(gridCIDs_sectionIDs_rootPixIDsRel,1); % including boundary cells
numGridCellFacesTot = numGridCellsTot * 6;

% Init
fm_cellsAll = readFmFile(pathToGridCellInteriorFeatures,1);
numGridCellFeatures = size(fm_cellsAll,2);
numCellFaceFeatures = numGridCellFeatures*3;
fm_cellFaces = zeros(numGridCellFacesTot,numCellFaceFeatures);
numGridCellsPerSlice = numGridCellsTot / numZ;
cellFaceIndStop_i = 0;
for i=1:numGridCellsTot
    if((sum(boundaryGridCellInds==i))==0)
        % get ID of direct neighbor for each face in the  given order
        cell_sectionID = floor(i/numGridCellsPerSlice);
        cellID_wrtSection = mod(i,numGridCellsPerSlice);
        [cellR,cellC] = ind2sub([numCellsY numCellsX],cellID_wrtSection);

        setOfDirectNeighbors_6 = getDirectFaceNeighborsInOrder...
                        (i,cell_sectionID,cellR,cellC,...
                        numCellsY,numCellsX,numZ);
        FMs_for_directNeighborSet = readFaceFeatures...
                        (fm_cellFaces,i,setOfDirectNeighbors_6);
        
        cellFaceIndStart_i = cellFaceIndStop_i +1;
        cellFaceIndStop_i = cellFaceIndStop_i +6;
        fm_cellFaces(cellFaceIndStart_i:cellFaceIndStop_i,:)...
                        = FMs_for_directNeighborSet;
    end
end

% save
fm_name = sprintf('fm_gridCellFacesAll.mat');
saveFileName = fullfile(saveFilePath,fm_name);
save(saveFileName,fm_cellFaces);