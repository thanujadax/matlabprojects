% create suggestions for the next iteration of 3DILP
function s_suggestions = createSuggestions()

% Inputs:
%   s_evaluation : output from the evaluation of the 3D segmentation
%   data : raw data

% Outputs:
%   new constraints - to avoid global constraint violations
%   new 2D hypothesis
%       - changes of priors ()
%           - for merge errors, change of state of voxels rewarding a
%           membrane through the 2D hypothesis with the merge error
%           - for split errors, change of state of voxels rewarding the
%           false membrane to be cell interior, combining the pair of cells
%           in concern to combine


[s_suggestions.A,s_suggestions.b,s_suggestions.senseArray]...
                    = getNewConstraints();


% getting 2D hypothesis : changes to the low-level costs

% 1 : merge errors
s_suggestions.membraneHypothesisGridCells ...
                    = getMembraneHypothesisGridCells();

% 2 : split errors
s_suggestions.cellInteriorHypothesisGridCells ... 
                    = getCellInteriorHypothesisGridCells();

