function seg3D = gridILP3D()
% 3D ILP on 3D grid
% Version 1.0
% 2014.02.25

%% Parameters
verbose = 2; % 0,1,2
usePrecomputedFeatureMatrices = 0;

gridResX = 4; % num pix
gridResY = 4;
gridResZ = 6; % distance between 2 adjacent slices in pixels

oriFiltLen = 29;    % window size for pix feature extraction and orientation filter
halfWidth_strucEl = 3; % structure_element_half_width-1 for orientation filter
csHist = oriFiltLen; % window size for histogram creation


numTrees = 500; % RFC

% weighting parameters for the cost function
W = [-10;    % cell interior
     -100;    % face 1
     -100;    % face 2
     -100;    % face 3
     -100;    % face 4
     -100;    % face 5
     -100;];  % face 6

NUM_VAR_PER_CELL = 7; 
% 1 - cell internal state
% 2 - face xy front
% 3 - face xy back
% 4 - face yz left
% 5 - face yz right
% 6 - face xz top
% 7 - face xz bottom



%% File paths
% 128x128 data set
pathForInputImages = '/home/thanuja/Dropbox/data/3D_Grid_ILP/stack1/raw/';
pathToFeatureMat = '/home/thanuja/Dropbox/data/3D_Grid_ILP/stack1/fm/';

% toy
% pathForInputImages = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/raw/';
% pathToFeatureMat = '/home/thanuja/Dropbox/data/3D_Grid_ILP/toy1/fm/';

% pathForInputImages ='/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/raw/';
% pathToFeatureMat ='/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/fm/';

pathToRFCs = '/home/thanuja/Dropbox/data/3D_Grid_ILP/trainingData/RFCs';

fileNameString = '*.png';
% subDir_Fm = 'fm'
subDir_sectionFm = 'indivSectionFMs';
subDir_cellInteriorFm = 'cellInteriorFMs';
subDir_cellFaceFm = 'cellFaceFMs';
subDir_RFCs = 'RFCs';

% fm file names
name_fm_cellInterior_sansBorder = 'fm_cellInterior_sansBorderCells.mat';
name_fm_cellInterior = 'fm_cellInterior.mat';
name_fm_faces12 = 'fm_faces12.mat';
name_fm_faces34 = 'fm_faces34.mat';
name_fm_faces56 = 'fm_faces56.mat';
% fm variable order
listInds_fm_cellsSansBorder_name = 'fm_listInds_cellsSansBorder.mat';
listInds_fm_cells_name = 'fm_listInds_cells.mat';
listInds_fm_face12_name = 'fm_listInds_face12.mat';
listInds_fm_face34_name = 'fm_listInds_face34.mat';
listInds_fm_face56_name = 'fm_listInds_face56.mat';

%% Read inputs
imageStack3D = readImages2StackWithBoundary...
                (pathForInputImages,fileNameString,gridResY,gridResX);
[numR,numC,numZ] = size(imageStack3D);

%% Create 3D grid
% addressing mechanism: gridID <--> pixels(of particular section)
% each gridCell correspond to a cube.
% each cube is made by expanding the corresponding anisotropic image patch
% vector: gridID|sectionID|pixels

% Membrane cubes are already added to the boundary
[gridCIDs_sectionIDs_rootPixIDsRel,numCellsY,numCellsX]...
            = getGridCells...
            (numR,numC,numZ,gridResX,gridResY,gridResZ);
%   griIDs_sectionIDs_rootPixIDsRel - matrix where each col is suggested by name
%   rootPixID is the pixInd of the start pixel (0,0,0) of each cell,
%   wrt each slice (relative coordinates in each slice)        
        
%   cellStats - number of grid cells along each dimension 
%       [numR,numC,numZ] = [numY,numX,numZ]
gridCellStats = [numCellsY, numCellsX,numZ];
disp('calculating border grid cell IDs')
[borderGridCellInds,borderFaceInds,borderCellFaceInds]...
                = getBoundaryCellInds(numCellsY,numCellsX,numZ);
disp('done')
numCellsPerSection = numCellsY * numCellsX;
numCells = numCellsY * numCellsX * numZ;
%% Compute feature matrices for gridCells and gridCellFaces
if(~usePrecomputedFeatureMatrices)
    % File names of feature matrices saved
%   <subDir_cellInteriorFm>/fm_cellInterior.mat
%   <subDir_cellFaceFm>/fm_faces12.mat      
%   <subDir_cellFaceFm>/fm_faces3456.mat  

% The face feature files exclude borderFaces and borderCells.

% calculate feature mats
% compute feature matrices for each section of the given stack of
% images, and save in the given path
disp('Calculating features for each section (pixelwise)')
computeFeaturesForEachSlice(pathToFeatureMat,subDir_sectionFm,imageStack3D,...
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
      
else
    % using precomputed feature matrices
    % check if available in the given path.
end
%% Unary activation scores from RFCs
unaryScoresMat = computeUnaryScoreMatRFC...
            (pathToFeatureMat,pathToRFCs,numTrees,numCells,NUM_VAR_PER_CELL,...
            subDir_cellInteriorFm,subDir_cellFaceFm,...
            name_fm_faces12,listInds_fm_face12_name,...
            name_fm_faces34,listInds_fm_face34_name,...
            name_fm_faces56,listInds_fm_face56_name,...
            name_fm_cellInterior,listInds_fm_cells_name);

%% ILP formulation
% constraints
[model.A,b,senseArray] = getILP3DGridConstraints(gridCellStats,borderGridCellInds);
% objective to minimized
f = getILP3DGridObjective(W,unaryScoresMat);

%% ILP solver
disp('using Gurobi ILP solver...');
% model.A = sparse(double(A));
model.rhs = b;
model.obj = f';
model.sense = senseArray;
% model.vtype = vtypeArray;
model.vtype = 'B';
% model.lb = lbArray;
% model.ub = ubArray;
model.modelname = '3D_Grid_ILP';
% initial guess
% model.start = labelVector;

params.LogFile = 'gurobi_3D_Grid_ILP.log';
params.Presolve = 0;
params.ResultFile = 'modelfile_3D_Grid_ILP.mps';
params.InfUnbdInfo = 1;

resultGurobi = gurobi(model,params);
x = resultGurobi.x;
%% Visualization
% numCellsPerSection = numCellsY * numCellsX;
rootPixels = gridCIDs_sectionIDs_rootPixIDsRel(1:numCellsPerSection,3);
simpleVisualizeStack(x,rootPixels,gridResY,gridResX,numR,numC,...
                    numCellsY,numCellsX,numZ);