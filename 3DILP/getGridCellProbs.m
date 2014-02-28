function gridCellInteriorUnary = getGridCellProbs...
            (fm,forestCellInteriorProb,numTrees)

% Inputs:
%   inputImageStackMat - 3D matrix containing adjacent slices of a stack,
%   including added boundary of cell exterior.
%   griIDs_sectionIDs_rootPixIDsRel - matrix where each col is suggested by name
%       rootPixID is the pixInd of the start pixel (0,0,0) of each cell,
%       wrt each slice (relative coordinates in each slice) 
%   forestCellInteriorProb
%   gridStats
%   numTrees

% Output:
%   gridCellInteriorUnary - the probabily that each gridCell is neuron
%   interior

% read precomputed featureMats for each image


% get features for gridCells

% Run classifier
[y_h,v] = classRF_predict(double(fm), forestCellInteriorProb);

gridCellInteriorUnary = v(:,2)./numTrees;