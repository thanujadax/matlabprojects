function createXYforRFtraining()
% create x y for RF training

% gridCellInterior labels - {0,1}
%   0 - interior
%   1 - exterior (membrane or outer)

% gridCellFace labels - {0,1}
%   0 - cell interior (same cell) on both sides or 
%           cell exterior on this side of the face (or both faces)
%   1 - cell interior of different cells on both sides or
%           next cell is exterior where this cell is interior

%% Parameters
verbose = 2; % 0,1,2
usePrecomputedFeatureMatrices = 0;

oriFiltLen = 29;    % window size for pix feature extraction and orientation filter
halfWidth_strucEl = 3; % structure_element_half_width-1 for orientation filter
csHist = oriFiltLen; % window size for histogram creation

doILPlabeling = 0;
fileNameString = '*.png';
gridResX = 128; % num pix
gridResY = 128;
gridResZ = 6; % distance between 2 adjacent slices in pixels

thresh_mem = 30; % membrane threshold % for initial labels for gridCells

NUM_VAR_PER_CELL = 7; 
% 1 - cell internal state
% 2 - face xy front
% 3 - face xy back
% 4 - face yz left
% 5 - face yz right
% 6 - face xz top
% 7 - face xz bottom

%% File paths and names
% % check if subdirectory exists. Create if not
% checkAndCreateSubDir(pathToFm,subDir_cellInteriorFm);

subDir_sectionFm = 'indivSectionFMs';
subDir_cellInteriorFm = 'cellInteriorFMs';

% rawImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/raw/';
% labelImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/neuron/';
% pathToFeatureMat = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/';
% saveLabelFilePath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/';

% toy
rawImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/raw/';
labelImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/labels/';
pathToFeatureMat = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/fm/';
saveLabelFilePath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/';

gridCellLabelsFileName = 'gridCellLabels.mat';
gridFace12LabelsFileName = 'gridCellFace12Labels.mat';
gridFace3456LabelsFileName = 'gridCellFace3456Labels.mat';

%% Read images
% raw
imageStack3D_raw = readImages2StackWithBoundary...
                (rawImgPath,fileNameString,gridResY,gridResX);
disp('Loaded raw images')
[sizeR,sizeC,numZ] = size(imageStack3D_raw);    

% neuron labels
imageStack3D_label = readImages2StackWithBoundary...
                (labelImgPath,fileNameString,gridResY,gridResX);
disp('Loaded training labels')
%% Create grid
% also, get initial labels for gridCellInteriors and gridCellFaces
disp('Creating grid with initial labels assigned to each grid cell..')
[gridCIDs_sectionIDs_rootPixIDsRel,gridCellInteriorInitLabels,...
    cellInteriorRGBLabelsAll,numCellsY,numCellsX]...
            = createInitLabelsForGridCells...
            (imageStack3D_label,sizeR,sizeC,numZ,gridResX,gridResY,gridResZ,...
            thresh_mem);
disp('done.')

disp('Assigning initial labels for gridCellFaces..')

gridCellFaceInitLabels = getCellFaceInitLabels...
                    (cellInteriorRGBLabelsAll,numZ,...
                    numCellsX,numCellsY);

gridCellStats = [numCellsY, numCellsX,numZ]; % used in ILP labeling
numCellsPerSection = numCellsY * numCellsX;
numCells = numCellsY * numCellsX * numZ;
% boundary cell ids
disp('calculating border grid cell IDs')
[borderGridCellInds,borderFaceInds]...
                = getBoundaryCellInds(numCellsY,numCellsX,numZ);

%% create x
% File names of feature matrices saved
%   <subDir_cellInteriorFm>/fm_cellInterior.mat
%   <subDir_cellFaceFm>/fm_faces12.mat      
%   <subDir_cellFaceFm>/fm_faces3456.mat  

% The face feature files exclude borderFaces and borderCells.

% calculate feature mats
% compute feature matrices for each section of the given stack of
% images, and save in the given path
disp('Calculating features for each section (pixelwise)')
computeFeaturesForEachSlice(pathToFeatureMat,subDir_sectionFm,imageStack3D_raw,...
            oriFiltLen, halfWidth_strucEl, csHist);
% writes one file (fm_gridCellInteriorAll.mat) with all gridCells of
% all sections
disp('Calculating features for grid cells ...')
computeFmGridCellInterior(pathToFeatureMat,subDir_cellInteriorFm,...
            subDir_sectionFm,numZ,gridCIDs_sectionIDs_rootPixIDsRel,...
            gridResX,gridResY);

disp('Calculating features for grid cell faces ...')        
computeFmGridFaces(pathToFeatureMat,borderGridCellInds,borderFaceInds,...
                gridCIDs_sectionIDs_rootPixIDsRel,numZ,numCellsY,numCellsX,...
                subDir_cellInteriorFm,subDir_cellFaceFm)
disp('********** Feature calculation done. Results saved. *********')

%% create y
% gridCellInteriorInitLabels
% gridCellFaceInitLabels

% visualize initial gridCell labels
disp('Visualizing initial grid cell labels..')
rootPixels = gridCIDs_sectionIDs_rootPixIDsRel(1:numCellsPerSection,3);
simpleVisualizeStack_activationVector...
                    (gridCellInteriorInitLabels,rootPixels,gridResY,gridResX,sizeR,sizeC,...
                    numCellsY,numCellsX,numZ);
% get labels from ILP
if(doILPlabeling)
    disp('ILP for structured label creation...')
    [gridCellInteriorLabels,gridCellFaceLabels,x]...
                        = getTrainingLabels3DgridILP...
                        (gridCellStats,borderGridCellIDs,...
                        gridCellInteriorInitLabels,gridCellFaceInitLabels);
end

%% save label mats
% save file 
% gridCellLabels sans border cells
saveFileName = fullfile(saveLabelFilePath,gridCellLabelsFileName);
if(doILPlabeling)
    gridCellInteriorLabels(borderGridCellInds) = [];
    save(saveFileName,'gridCellInteriorLabels');
else    
    gridCellInteriorInitLabels(borderGridCellInds) = [];
    save(saveFileName,'gridCellInteriorInitLabels');
end
disp('Saved labels for grid cells at:')
disp(saveFileName)

% gridCellFaceLabels sans border faces
% convert into vector
numFaces = numCells * 6;
faceLabels = reshape(gridCellFaceInitLabels',numFaces,1);

% 1. faces12
numXYfaces = numCells * 2;
seq = 1:numCells;
face1IDs = (seq-1)*6 + 1;
face2IDs = (seq-1)*6 + 2;
face12IDs = [face1IDs face2IDs];
% remove border face ids
face12IDs = setdiff(face12IDs,borderFaceInds);
labels_faces12 = faceLabels(face12IDs);
saveFileName = fullfile(saveLabelFilePath,gridFace12LabelsFileName);
save(saveFileName,'labels_faces12');
disp('Saved labels for grid cell faces12 at:')
disp(saveFileName)

% 2. faces3456
seq2 = 1: numCells*6;
face3456IDs = setdiff(seq2,face12IDs);
face3456IDs = setdiff(face3456IDs,borderFaceInds);
labels_faces3456 = faceLabels(face3456IDs);

saveFileName = fullfile(saveLabelFilePath,gridFace3456LabelsFileName);
save(saveFileName,'labels_faces3456');
disp('Saved labels for grid cell faces3456 at:')
disp(saveFileName)

