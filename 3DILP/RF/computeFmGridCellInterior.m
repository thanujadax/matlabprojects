function computeFmGridCellInterior(pathToFm,subDir_cellInteriorFm,...
                subDir_sectionFm,numZ,gridCIDs_sectionIDs_rootPixIDsRel,...
                gridResX,gridResY)

% computes and saves the feature matrices for each gridCell

% Inputs:
%   numZ - number of sections in the stack


% path to save section features
saveFilePath = fullfile(pathToFm,subDir_cellInteriorFm);
pathToSliceFeatures = fullfile(pathToFm,subDir_sectionFm);

% check if subdirectory exists. Create if not
checkAndCreateSubDir(pathToFm,subDir_cellInteriorFm);

% featureFile.mat name structure: fm_slice_%d.mat


numGridCellsTot = size(gridCIDs_sectionIDs_rootPixIDsRel,1); % including boundary cells

% Init
onlyFile = 0;
fileIndex = 1;
fm_slice = readFmFile(pathToSliceFeatures,fileIndex);
[sizeR,sizeC,numPixFeatures] = size(fm_slice);
numGridCellFeatures = numPixFeatures * 5; % min,max,mean,var,medain
fm_cellInterior = zeros(numGridCellsTot,numGridCellFeatures);
numGridCellsPerSlice = numGridCellsTot / numZ;
% gridIndStop_slice_i = 0;

for i=1:numZ
    % read slice features
    fileIndex = i;
    fm_slice = readFmFile(pathToSliceFeatures,fileIndex);
    % get the features relevant for each gridCell in the precomputed order
    % get the rootGridCellPix for the current section
    tmp_rowInd_logical = (gridCIDs_sectionIDs_rootPixIDsRel(:,2)==i);
    rootPixIndsOfCells = gridCIDs_sectionIDs_rootPixIDsRel(tmp_rowInd_logical,3);
    gridCellInds = gridCIDs_sectionIDs_rootPixIDsRel(tmp_rowInd_logical,1);
    gridCellInds = gridCellInds + (i-1) * numGridCellsPerSlice;
    fm_gridsForSlice_i = convertPixFmToCellFm...
        (fm_slice,rootPixIndsOfCells,gridResX,gridResY,sizeR,sizeC,...
        numGridCellsPerSlice);
    % add to feature matrix fm_cellInterior
    % gridIndStart_slice_i = (i-1) * numGridCellsPerSection +1;
    % gridIndStop_slice_i = gridIndStart_slice_i + numGridCellsPerSection -1;
    fm_cellInterior(gridCellInds,:) = fm_gridsForSlice_i;
    clear fm_gridsForSlice_i
end

% save
fm_name = sprintf('fm_gridCellInteriorAll.mat');
saveFileName = fullfile(saveFilePath,fm_name);
save(saveFileName,fm_cellInterior);
