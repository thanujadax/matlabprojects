function s_3DcompEval = compEval3D(neuronComp3Dfeatures,s_sketetonEval)

% Inputs:
%   neuronComp3D
%   skeletonEvaluation
%       loop links : global constraint vilations
%       abrupt ends
%       locally confident suggraphs
%       locally weak subgraphs

% Outputs:
%   link strength
%       links with very low strength should be turned off

% calculate features

% get potential 3D links for abrupt ends

features_3Dcomps = getFeatures_3Dcomps();
s_3DcompEval.linkScores = RFC_getLinkScores(features_3Dcomps);

