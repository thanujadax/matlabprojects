function s_3DcompEval = compEval3D(neuronComp3Dfeatures,s_sketetonEval)

% Inputs:
%   neuronComp3D
%   skeletonEvaluation

% Outputs:
%   weak 3D links
%   strong 3D links

% calculate features

features_3Dcomps = getFeatures_3Dcomps();
s_3DcompEval.linkScores = RFC_getLinkScores(features_3Dcomps);
