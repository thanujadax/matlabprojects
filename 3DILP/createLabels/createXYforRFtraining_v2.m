function createXYforRFtraining_v2()
% create x y for RF training

% gridCellInterior labels - {0,1}
%   0 - interior
%   1 - exterior (membrane + mitochondria)***

% gridCellFace labels - {0,1}
%   0 - cell interior (same cell) on both sides or 
%           cell exterior on this side of the face (or both faces)
%   1 - cell interior of different cells on both sides or
%           next cell is exterior where this cell is interior

%% Parameters
oriFiltLen = 29;    % window size for pix feature extraction and orientation filter
halfWidth_strucEl = 3; % structure_element_half_width-1 for orientation filter
csHist = oriFiltLen; % window size for histogram creation

doILPlabeling = 0;
fileNameString = '*.png';
gridResX = 4; % num pix
gridResY = 4;
gridResZ = 6; % distance between 2 adjacent slices in pixels

thresh_mem = 40; % membrane threshold % for initial labels for gridCells

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
subDir_cellFaceFm = 'cellFaceFMs';
subDir_RFCs = 'RFCs';

rawImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/testData3/raw/';
neuronLabelImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/testData3/neuron/';
mitoLabelImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/testData3/mito/';

pathToFeatureMat = '/home/thanuja/Dropbox/data/3D_Grid_ILP/testData3/fm/';
saveLabelFilePath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/testData3/lblmat';

% toy
% rawImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/raw/';
% labelImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/labels/';
% pathToFeatureMat = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/fm2/';
% saveLabelFilePath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/lblmat2';

gridCellLabelsFileName = 'gridCellLabels.mat';
gridFace12LabelsFileName = 'gridCellFace12Labels.mat';
gridFace34LabelsFileName = 'gridCellFace34Labels.mat';
gridFace56LabelsFileName = 'gridCellFace56Labels.mat';

listInds_fm_cellsSansBorder_name = 'fm_listInds_cellsSansBorder.mat';
listInds_fm_cells_name = 'fm_listInds_cells.mat';
listInds_fm_face12_name = 'fm_listInds_face12.mat';
listInds_fm_face34_name = 'fm_listInds_face34.mat';
listInds_fm_face56_name = 'fm_listInds_face56.mat';

% fm file names
name_fm_cellInterior_sansBorder = 'fm_cellInterior_sansBorderCells.mat';
name_fm_cellInterior = 'fm_cellInterior.mat';
name_fm_faces12 = 'fm_faces12.mat';
name_fm_faces34 = 'fm_faces34.mat';
name_fm_faces56 = 'fm_faces56.mat';

%% Read images
% raw
imageStack3D_raw = readImages2StackWithBoundary...
                (rawImgPath,fileNameString,gridResY,gridResX);
disp('Loaded raw images')
[sizeR,sizeC,numZ] = size(imageStack3D_raw);    

% neuron labels
imageStack3D_label_neuron = readImages2StackWithBoundary...
                (neuronLabelImgPath,fileNameString,gridResY,gridResX);
            
% mito labels
imageStack3D_label_mito = readImages2StackWithBoundary...
                (mitoLabelImgPath,fileNameString,gridResY,gridResX);
disp('Loaded training labels')
%% Create grid
% also, get initial labels for gridCellInteriors and gridCellFaces
disp('Creating grid with initial labels assigned to each grid cell..')
[gridCIDs_sectionIDs_rootPixIDsRel,gridCellInteriorInitLabels,...
    cellInteriorRGBLabelsAll,numCellsY,numCellsX]...
            = createInitLabelsForGridCells...
            (imageStack3D_label_neuron,sizeR,sizeC,numZ,gridResX,gridResY,gridResZ,...
            thresh_mem,imageStack3D_label_mito);
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
[borderGridCellInds,borderFaceInds,borderCellFaceInds]...
                = getBoundaryCellInds(numCellsY,numCellsX,numZ);
disp('done')
%% Get feature matrices
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
            gridResX,gridResY,borderGridCellInds,listInds_fm_cells_name,...
            listInds_fm_cellsSansBorder_name,name_fm_cellInterior,...
                name_fm_cellInterior_sansBorder);

disp('Calculating features for grid cell faces ...')        
computeFmGridFaces(pathToFeatureMat,borderGridCellInds,borderCellFaceInds,...
                gridCIDs_sectionIDs_rootPixIDsRel,numZ,numCellsY,numCellsX,...
                subDir_cellInteriorFm,subDir_cellFaceFm,...
    listInds_fm_face12_name,listInds_fm_face34_name,listInds_fm_face56_name,...
    name_fm_faces12,name_fm_faces34,name_fm_faces56,name_fm_cellInterior);
disp('********** Feature calculation done. Results saved. *********')
% fm_faces12.mat
% fm_faces34.mat
% fm_faces56.mat
%% Get label vectors
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
faceLIDdir = fullfile(pathToFeatureMat,subDir_cellFaceFm);

% 1. faces12
face12LIDfile = fullfile(faceLIDdir,listInds_fm_face12_name);
face12LIDs = importdata(face12LIDfile);
labels_faces12 = faceLabels(face12LIDs);

saveFileName = fullfile(saveLabelFilePath,gridFace12LabelsFileName);
save(saveFileName,'labels_faces12');
disp('Saved labels for grid cell faces12 at:')
disp(saveFileName)

% 2. faces xz {3,4,3,4,...}
face34LIDfile = fullfile(faceLIDdir,listInds_fm_face34_name);
face34LIDs = importdata(face34LIDfile);
labels_faces34 = faceLabels(face34LIDs);

saveFileName = fullfile(saveLabelFilePath,gridFace34LabelsFileName);
save(saveFileName,'labels_faces34');
disp('Saved labels for grid cell faces34 at:')
disp(saveFileName)


% 3. faces yz {5,6,5,6,...}
face56LIDfile = fullfile(faceLIDdir,listInds_fm_face56_name);
face56LIDs = importdata(face56LIDfile);
labels_faces56 = faceLabels(face56LIDs);

saveFileName = fullfile(saveLabelFilePath,gridFace56LabelsFileName);
save(saveFileName,'labels_faces56');
disp('Saved labels for grid cell faces56 at:')
disp(saveFileName)

disp('*******label vectors saved *******')

