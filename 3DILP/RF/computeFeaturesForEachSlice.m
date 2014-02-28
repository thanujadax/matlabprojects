function computeFeaturesForEachSlice(pathToFm,subDir_sectionFm,imageStack3D,...
                oriFiltLen, halfWidth_strucEl, csHist)

% computes and saves the feature matrices for each section of the input
% image stack. fm(i,j,k) corresponds to the feature k of the gridPatch
% (i,j) of each section.

% path to save section features
saveFilePath = fullfile(pathToFm,subDir_sectionFm);

% check if subdirectory exists. Create if not
checkAndCreateSubDir(pathToFm,subDir_sectionFm);

% featureFile.mat name structure: fm_slice_%d.mat

[~,~,numZ] = size(imageStack3D);

for i=1:numZ
    % get feature matrix for this section
    fm  = getFeatures_nonMembrane...
            (imageStack3D(:,:,i), oriFiltLen, halfWidth_strucEl, csHist);
    % save 
    fm_name = sprintf('fm_slice_%d.mat',i);
    saveFile = fullfile(saveFilePath,fm_name);
    save(saveFile,fm);
end