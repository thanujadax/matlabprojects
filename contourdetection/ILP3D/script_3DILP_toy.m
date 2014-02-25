% 3D ILP toy problem
% 2014.02.19

% ver 1.0
% Uses undirected edges on each 2D slice

%% Param
% 2D ILP
w_off_e = 16.8913;
w_on_e = -16.8913;
w_off_n = -10.6301;
w_on_n = 9.27912;
w_off_r = 1.07527;
w_on_r = -1.07527;

%%  Toy inputs (2D grids)
%   with appropriate feature values for edges, regions and triangles
%   2D variables (edges,nodes,regions) for each section are stored in a
%   cell array
%   3D variables (triangles) are organized in a cell array, each cell
%   corresponding to a pair of adjacent sections

%% read inputs sections


%% extract the data structures for ILP


%%  ILP formulation

c_TsetsWithCommonEdges = [];  
% not required as long as lateral contiguity is not a constraint
            
[model.A,b,senseArray,numEdges,numNodeConf,numRegions,nodeTypeStats] = ...
            getILPConstraints3Dtoy...
            (c_edgeIDs,c_Ts_on,sectVarStats,c_TsetsWithCommonEdges,...
            edgeListInds,edges2nodes,nodeEdgeIDs,junctionTypeListInds,...
            jEdges,dirEdges2regionsOnOff,setOfRegions,activeWSregionListInds_tr);
        
f = getILPObj3Dtoy(edgeUnary,nodeAngleCosts,...
            regionUnary,T_unary,...
            w_off_e,w_on_e,w_off_n,w_on_n,w_off_r,w_on_r,w_on_T,w_off_T);

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
