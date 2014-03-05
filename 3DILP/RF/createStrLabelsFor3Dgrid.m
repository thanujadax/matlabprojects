% function createStrLabelsFor3Dgrid(rawImgPath,labelImgPath,saveLabelFilePath,...
%                 fileNameString,gridResX,gridResY,gridResZ)

% create structured labels for

% Inputs:
%   rawImgPath
%   labelImgPath
%   saveLabelFilePath

% gridCellInterior labels - {0,1}
%   0 - interior
%   1 - exterior (membrane or outer)

% gridCellFace labels - {0,1}
%   0 - cell interior (same cell) on both sides or 
%           cell exterior on this side of the face (or both faces)
%   1 - cell interior of different cells on both sides or
%           next cell is exterior where this cell is interior

disp('Creating structured training labels for 3DgridILP')
%% File paths
rawImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/raw/';
labelImgPath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/neuron/';
saveLabelFilePath = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/';
gridCellLabelsFileName = 'gridCellLabels.mat';
gridCellFaceLabelsFileName = 'gridCellFaceLabels.mat';

%% Params
fileNameString = '*.png';
gridResX = 16; % num pix
gridResY = 16;
gridResZ = 6; % distance between 2 adjacent slices in pixels

thresh_mem = 30; %

%% Read inputs
% Load images into 3D matrix
% Load labels into 3D matrix
% add boundary cells to both stacks
imageStack3D_raw = readImages2StackWithBoundary...
                (rawImgPath,fileNameString,gridResY,gridResX);
disp('Loaded raw images')
[sizeR,sizeC,numZ] = size(imageStack3D_raw);
% neuron labels
imageStack3D_label = readImages2StackWithBoundary...
                (labelImgPath,fileNameString,gridResY,gridResX);
disp('Loaded training labels')
            
%% Create grid
% also, get initial labels for gridCellInteriors
disp('Creating grid with initial labels assigned to each grid cell..')
[gridIDs_sectionIDs_rootPixIDsRel,gridCellInteriorInitLabels,...
            numCellsY,numCellsX]...
            = createInitLabelsForGridCells...
            (imageStack3D_label,sizeR,sizeC,numZ,gridResX,gridResY,gridResZ,...
            thresh_mem);
disp('done.')

gridCellStats = [numCellsY, numCellsX,numZ];


%% Get cell interior labels for all cells
% boundary cell ids
disp('calculating border grid cell IDs')
borderGridCellIDs = getBoundaryCellInds(numCellsY,numCellsX,numZ);
% visualize initial gridCell labels
disp('Visualizing initial grid cell labels..')
rootPixels = gridIDs_sectionIDs_rootPixIDsRel(1:numCellsPerSection,3);
simpleVisualizeStack_activationVector...
                    (gridCellInteriorInitLabels,rootPixels,gridResY,gridResX,sizeR,sizeC,...
                    numCellsY,numCellsX,numZ);
% get labels from ILP
disp('ILP for structured label creation...')
[gridCellInteriorLabels,gridCellFaceLabels,x]...
                    = getTrainingLabels3DgridILP...
                    (gridCellStats,borderGridCellIDs,gridCellInteriorInitLabels);
%% Visualize
numCellsPerSection = numCellsY * numCellsX;

simpleVisualizeStack(x,rootPixels,gridResY,gridResX,sizeR,sizeC,...
                    numCellsY,numCellsX,numZ)
%% Get face labels for all 6 faces of each cell
% save file 
% gridCellLabels
saveFileName = fullfile(saveLabelFilePath,gridCellLabelsFileName);
save(saveFileName,gridCellInteriorLabels);
disp('Saved labels for grid cells at:')
disp(saveFileName)
% gridCellFaceLabels
saveFileName = fullfile(saveLabelFilePath,gridCellFaceLabelsFileName);
save(saveFileName,gridCellFaceLabels);
disp('Saved labels for grid cell faces at:')
disp(saveFileName)