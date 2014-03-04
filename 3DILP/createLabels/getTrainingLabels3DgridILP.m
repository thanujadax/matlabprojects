function [gridCellInteriorLabels,gridCellFaceLabels]...
                    = getTrainingLabels3DgridILP...
                    (cellStats,borderCellIDs,gridCellInteriorInitLabels)

% Inputs:
%   cellStats - number of grid cells along each dimension 
%       [numR,numC,numZ] = [numY,numX,numZ]


% Outputs:




%% ILP formulation
disp('Formulating constraints...')
% constraints
numBorderCells = numel(borderCellIDs);
[model.A,b,senseArray] = getILP3DGridConstraints(cellStats,numBorderCells);
disp('done.')
% objective to minimize
f = getILP3DGridObjective_traindingLabels(gridCellInteriorInitLabels,borderCellIDs);

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
model.modelname = '3D_Grid_ILP_trainigLabels';
% initial guess
% model.start = labelVector;

params.LogFile = 'gurobi_3D_Grid_ILP.log';
params.Presolve = 0;
params.ResultFile = 'modelfile_3D_Grid_ILP.mps';
params.InfUnbdInfo = 1;

resultGurobi = gurobi(model,params);
x = resultGurobi.x;
%% Extract labels separately: gridCellInterior and gridCellFaces
% Column structure of A (structure of variable vector
%   Each grid cell gives rise to 7 state variables
%   1 - cell internal state
%   2 - state of face x-y (front)
%   3 - state of face x-y (back)
%   4 - state of face y-z (left)
%   5 - state of face y-z (right)
%   6 - state of face x-z (top)
%   7 - state of face x-z (bottom)
numCells = numel(gridCellInteriorInitLabels);
outLabels = reshape(x,7,numCells);
% take transpose
outLabels = outLabels';
gridCellInteriorLabels = outLabels(:,1);
gridCellFaceLabels = outLabels(:,2:7);

