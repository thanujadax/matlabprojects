function computeFmGridCellInterior(pathToFm,subDir_cellInteriorFm,...
                subDir_sectionFm,imageStack3D,oriFiltLen,...
                halfWidth_strucEl, csHist)

% computes and saves the feature matrices for each gridCell


% path to save section features
saveFilePath = fullfile(pathToFm,subDir_cellInteriorFm);
pathToSliceFeatures = fullfile(pathToFm,subDir_sectionFm);

% check if subdirectory exists. Create if not
checkAndCreateSubDir(pathToFm,subDir_cellInteriorFm);

% featureFile.mat name structure: fm_slice_%d.mat

[~,~,numZ] = size(imageStack3D);

% Init
fm_cellInterior = zeros(numGridCells,numFeatures);

for i=1:numZ
    % read slice features
    fm_slice = readSliceFeatures(i);
    % get the features relevant for each gridCell in the precomputed order
    fm_gridsForSlice_i = convertPixFmToCellFm(fm_slice);
    % add to feature matrix fm_cellInterior
    gridIndStart
    gridIndStop
    fm_cellInterior(gridIndStart:gridIndStop,:) = fm_gridsForSlice_i;
    clear fm_gridsForSlice_i
    
end

% save
fm_name = sprintf('fm_gridCellInteriorAll.mat');
saveFile = fullfile(saveFilePath,fm_name);
save(saveFile,fm_cellInterior);
