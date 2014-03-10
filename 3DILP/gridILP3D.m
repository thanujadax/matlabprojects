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

NUM_VAR_PER_CELL = 7; 
% 1 - cell internal state
% 2 - face xy front
% 3 - face xy back
% 4 - face yz left
% 5 - face yz right
% 6 - face xz top
% 7 - face xz bottom

%% File paths
pathForInputImages = '/home/thanuja/Dropbox/data/3D_Grid_ILP/stack1/';
fileNameString = '*.png';

pathToFeatureMat = '/home/thanuja/Dropbox/data/3D_Grid_ILP/stack1/fm/';
subDir_sectionFm = 'indivSectionFMs';
subDir_cellInteriorFm = 'cellInteriorFMs';
subDir_cellFaceFm = 'cellFaceFMs';
subDir_cellFaceFMs_xy1 = 'cellFaceFMs_xy1';
subDir_cellFaceFMs_xy2 = 'cellFaceFMs_xy2';
subDir_cellFaceFMs_yz1 = 'cellFaceFMs_yz1';
subDir_cellFaceFMs_yz2 = 'cellFaceFMs_yz2';
subDir_cellFaceFMs_xz1 = 'cellFaceFMs_xz1';
subDir_cellFaceFMs_xz2 = 'cellFaceFMs_xz2';

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
[boundaryGridCellInds,borderFaceInds] ...
        = getBoundaryCellInds(numCellsY,numCellsX,numZ);

%% Compute feature matrices for gridCells and gridCellFaces
if(~usePrecomputedFeatureMatrices)
    % compute feature matrices for each section of the given stack of
    % images, and save in the given path
    computeFeaturesForEachSlice(pathToFeatureMat,subDir_sectionFm,imageStack3D,...
                oriFiltLen, halfWidth_strucEl, csHist);
    
    % writes one file (fm_gridCellInteriorAll.mat) with all gridCells of
    % all sections
    computeFmGridCellInterior(pathToFeatureMat,subDir_cellInteriorFm,...
                subDir_sectionFm,numZ,gridCIDs_sectionIDs_rootPixIDsRel,...
                gridResX,gridResY);
    
    computeFmGridFaces(pathToFeatureMat,boundaryGridCellInds,...
                gridCIDs_sectionIDs_rootPixIDsRel,numZ,numCellsY,numCellsX);
      
else
    % using precomputed feature matrices
    % check if available in the given path.
end
%% Unary activation scores from RFCs
unaryScoresMat = computeUnaryScoreMatRFC();

%% ILP formulation

% constraints
[model.A,b,senseArray] = getILP3DGridConstraints(gridCellStats);
% objective to minimize
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
seg3D = [];