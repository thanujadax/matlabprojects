function computeFmGridCellInterior(pathToFm,subDir_cellInteriorFm,...
                subDir_sectionFm,imageStack3D,oriFiltLen,...
                halfWidth_strucEl, csHist,gridCIDs_sectionIDs_rootPixIDsRel)

% computes and saves the feature matrices for each gridCell


% path to save section features
saveFilePath = fullfile(pathToFm,subDir_cellInteriorFm);
pathToSliceFeatures = fullfile(pathToFm,subDir_sectionFm);

% check if subdirectory exists. Create if not
checkAndCreateSubDir(pathToFm,subDir_cellInteriorFm);

% featureFile.mat name structure: fm_slice_%d.mat

[~,~,numZ] = size(imageStack3D);
numGridCellsTot = size(gridCIDs_sectionIDs_rootPixIDsRel,1); % including boundary cells

% Init
fm_cellInterior = zeros(numGridCellsTot,numFeatures);
% gridIndStop_slice_i = 0;
for i=1:numZ
    % read slice features
    fm_slice = readSliceFeatures(pathToSliceFeatures,i);
    % get the features relevant for each gridCell in the precomputed order
    fm_gridsForSlice_i = convertPixFmToCellFm(fm_slice);
    % add to feature matrix fm_cellInterior
    gridIndStart_slice_i = (i-1) * numGridCellsPerSection +1;
    gridIndStop_slice_i = gridIndStart_slice_i + numGridCellsPerSection -1;
    fm_cellInterior(gridIndStart_slice_i:gridIndStop_slice_i,:)...
                = fm_gridsForSlice_i;
    clear fm_gridsForSlice_i
    
end

% save
fm_name = sprintf('fm_gridCellInteriorAll.mat');
saveFileName = fullfile(saveFilePath,fm_name);
save(saveFileName,fm_cellInterior);
