% 3D ILP toy problem

%%  Define set of toy slices (2D grids)
%   with appropriate feature values for edges, regions and triangles
%   2D variables (edges,nodes,regions) for each section are stored in a
%   cell array
%   3D variables (triangles) are organized in a cell array, each cell
%   corresponding to a pair of adjacent sections


%%  ILP formulation
[model.A,b,senseArray] = getILPConstraints3Dtoy();
f = getILPObj3Dtoy();

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
model.modelname = 'ILP3Dtoy_1';
% initial guess
% model.start = labelVector;

params.LogFile = 'gurobi.log';
params.Presolve = 0;
params.ResultFile = 'modelfile.mps';
params.InfUnbdInfo = 1;

resultGurobi = gurobi(model,params);
x = resultGurobi.x;
%% visualization
