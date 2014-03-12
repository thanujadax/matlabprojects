function computeFmGridCellInterior(pathToFm,subDir_cellInteriorFm,...
                subDir_sectionFm,numZ,gridCIDs_sectionIDs_rootPixIDsRel,...
                gridResX,gridResY,borderCellIDs,listInds_fm_cells_name,...
                listInds_fm_cellsSansBorder_name,name_fm_cellInterior,...
                name_fm_cellInterior_sansBorder)

% computes and saves the feature matrices for each gridCell
% excludes border cells

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
numGridCellFeatures = numPixFeatures * 5; % min,max,mean,var,medain of each
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
    % gridCellInds = gridCellInds + (i-1) * numGridCellsPerSlice;
    fm_gridsForSlice_i = convertPixFmToCellFm...
        (fm_slice,rootPixIndsOfCells,gridResX,gridResY,sizeR,sizeC,...
        numGridCellsPerSlice,numGridCellFeatures);
    nanInds = isnan(fm_gridsForSlice_i);
    if(sum(sum(nanInds))>0)
        fm_gridsForSlice_i(nanInds) = 0;
    end
    % add to feature matrix fm_cellInterior
    % gridIndStart_slice_i = (i-1) * numGridCellsPerSection +1;
    % gridIndStop_slice_i = gridIndStart_slice_i + numGridCellsPerSection -1;
    fm_cellInterior(gridCellInds,:) = fm_gridsForSlice_i;
    clear fm_gridsForSlice_i
end

% save
% fm_name = sprintf('fm_cellInterior.mat');
saveFileName = fullfile(saveFilePath,name_fm_cellInterior);
save(saveFileName,'fm_cellInterior');

% save another version without the border cells
% remove the rows corresponding to borderCellIDs
fm_cellInterior_sansBorderCells = fm_cellInterior;
fm_cellInterior_sansBorderCells(borderCellIDs,:) = [];
% fm_name = sprintf('fm_cellInterior_sansBorderCells.mat');
saveFileName = fullfile(saveFilePath,name_fm_cellInterior_sansBorder);
save(saveFileName,'fm_cellInterior_sansBorderCells');

% save the cellListInds of which the features are calculated, in the same
% order as in the fm file
seq = 1:numGridCellsTot;
saveFileName = fullfile(saveFilePath,listInds_fm_cells_name);
save(saveFileName,'seq');
% sans border cells
listInds_fm_cells_sansBorder = setdiff(seq,borderCellIDs);
saveFileName = fullfile(saveFilePath,listInds_fm_cellsSansBorder_name);
save(saveFileName,'listInds_fm_cells_sansBorder');

