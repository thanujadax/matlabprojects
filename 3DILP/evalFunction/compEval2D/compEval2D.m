function s_2DcompEval = compEval2D(neuronComp2D,s_3DcompEval,s_skeletonEval)

% Inputs
%   neuronComp2D
%   s_skeletonEval
%       loop links : global constraint vilations
%       abrupt ends
%       locally confident suggraphs
%       locally weak subgraphs
%   s_3DcompEval
%       weak 3D links
%       strong 3D links

% Outputs 
%   potential merge errors (node IDs)
%   potential split errors (sets of nodeIDs)
%   potential missclassification (cellInterior -> membrane or vice versa)


comp2DFeatures = get2DcompFeatures();
s_2DcompEval.comp2Dscores = RFC_2Dcomponents(comp2DFeatures);

s_2DcompEval.mergeErrors = RFC_mergeErrors();
s_2DcompEval.splitErrors = RFC_splitErrors();

