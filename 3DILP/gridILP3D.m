function seg3D = gridILP3D()
% 3D ILP on 3D grid
% Version 1.0
% 2014.02.25

%% Parameters
verbose = 2; % 0,1,2
pathForInputImages = '/home/thanuja/Dropbox/data/3D_Grid_ILP/stack1/';
fileNameString = '*.png';
gridResX = 4; % num pix
gridResY = 4;
gridResZ = 6; % distance between 2 adjacent slices in pixels

NUM_VAR_PER_CELL = 7; 
% 1 - cell internal state
% 2 - face xy front
% 3 - face xy back
% 4 - face yz left
% 5 - face yz right
% 6 - face xz top
% 7 - face xz bottom

%% Read inputs
inputSections = dir(strcat(pathForInputImages,fileNameString));
numInputSections = length(inputSections);
% load images to 3D matrix
% read the first file to read the dimensions
inputImage1_FilePath = fullfile(pathForInputImages,inputSections(1).name);
im1 = uint8(imread(inputImage1_FilePath)); % uint8,single,double are the available options
[numR,numC] = size(im1);
clear im1
% Adjusting image dimensions to include boundary (membrane) cells
numR = numR + gridResY*2;
numC = numC + gridResX*2;
numZ = numInputSections + 2; % Thickness of slices not yet considered
inputImageStackMat = uint8(ones(numR,numC,numZ));
rStart = gridResY + 1;
rStop = numR - gridResY;
cStart = gridResX + 1;
cStop = numC - gridResX;
for i=2:numInputSections
    % section 1 is already filled with 1s (membrane)
    inputImageFilePath = fullfile(pathForInputImages,inputSections(i).name);
    inputImageStackMat(rStart:rStop,cStart:cStop,i) = ...
                uint8(imread(inputImageFilePath)); 
    % uint8,single,double are the available options
end

%% Create 3D grid
% addressing mechanism: gridID <--> pixels(of particular section)
% each gridCell correspond to a cube.
% each cube is made by expanding the corresponding anisotropic image patch
% vector: gridID|sectionID|pixels

% Membrane cubes are added to the boundary
griIDs_sectionIDs_rootPixIDsRel = getGridCells...
            (numR,numC,numZ,gridResX,gridResY,gridResZ);
%   griIDs_sectionIDs_rootPixIDsRel - matrix where each col is suggested by name
%   rootPixID is the pixInd of the start pixel (0,0,0) of each cell,
%   wrt each slice (relative coordinates in each slice)        
        
numCells = size(griIDs_sectionIDs_rootPixIDsRel,1);
%% Unary activation scores from RFCs
unaryScoreMat = zeros(numCells,NUM_VAR_PER_CELL);

unaryScoreMat(:,1) = getGridCellProbs(); % gridCellProbs_interior
unaryScoreMat(:,2) = getGridCellFaceProbs_xy(); % gridCellFaceProbs_xy_1 
unaryScoreMat(:,3) = getGridCellFaceProbs_xy(); % gridCellFaceProbs_xy_2 
unaryScoreMat(:,4) = getGridCellFaceProbs_sides(); % gridCellFaceProbs_yz_3 
unaryScoreMat(:,5) = getGridCellFaceProbs_sides(); % gridCellFaceProbs_yz_4 
unaryScoreMat(:,6) = getGridCellFaceProbs_sides(); % gridCellFaceProbs_xz_5 
unaryScoreMat(:,7) = getGridCellFaceProbs_sides(); % gridCellFaceProbs_xz_6 

%% ILP formulation

% constraints
[model.A,b,senseArray] = getILP3DGridConstraints(cellStats);
% objective to minimize
f = getILP3DGridObjective(W,unaryScoreMat);

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