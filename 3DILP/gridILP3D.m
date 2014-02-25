function seg3D = gridILP3D
% 3D ILP on 3D grid
% Version 1.0
% 2014.02.25

%% Parameters
verbose = 2; % 0,1,2
pathForInputImages = '';
fileNameString = '*.tif';


%% Read inputs
inputSections = dir(strcat(pathForImages_training,fileNameString)); 
%% Create 3D grid

%% Unary activation scores from RFCs

%% ILP formulation

% constraints
[model.A,b,senseArray] = getILP3DGridConstraints(cellStats);
% objective to minimize
f = getILP3DGridObjective(W,cellUnaries);

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
