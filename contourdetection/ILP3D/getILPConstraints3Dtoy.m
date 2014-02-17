function [A,b,senseArray] = getILPConstraints3Dtoy(c_edgeIDs)

% 2014.02.17

% Outputs:
%   A - sparse matrix containing coefficients for ILP
%   b - 
%   senseArray -

% Inputs

%% Init
rowStop = 0;
numSections = size(c_edgeIDs,2);

totRowsA = 
totColsA = 
A = sparse(totRowsA,totColsA);
%% 2D constraints
% formulated for each of the sections

%% 3D constraints
% formulated for each adjacent pair of sections


