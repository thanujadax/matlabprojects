function s_2DcompEval = compEval2D(neuronComp2D,s_3DcompEval,s_skeletonEval)

% Inputs
%   neuronComp2D
%   s_3DcompEval
%   s_skeletonEval
%       global constraint vilations
%       abrupt ends

% Outputs 
%   potential merge errors (node IDs)
%   potential split errors (sets of nodeIDs)


comp2Dfeatures = get2DcompFeatures();
s_2DcompEval.comp2Dscores = RFC_2Dcomponents(comp2DFeatures);

s_2DcompEval.mergeErrors = RFC_mergeErrors();
s_2DcompEval.splitErrors = RFC_splitErrors();

