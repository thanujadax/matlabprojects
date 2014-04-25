function s_skeletonEval = skeletonEval(adjGraph)

% Evaluates the skeleton structure

% Inputs:
%   precalculated skeleton features for classifier


% Outputs:
%   s_skeletonEval.loopLinks : global constraint violations
%       - mutually exclusive sets of links (linkIDs)

%   s_skeletonEval.twoWayMerges
%   s_skeletonEval.twoWaySplits

%   s_skeletonEval.abruptEnds (nodeIDs)

%   s_skeletonEval.locallyConfidentSubgraphs : (linkIDs)
%       - medium range subgraphs
%   s_skeletonEval.locallyWeakSubgraphs : (linkIDS)
%       - weak medium range subgraphs (linkIDs)

% global constraint violations
[s_skeletonEval.loopLinks, ...
  s_skeletonEval.twoWayMerges, ...
  s_skeletonEval.twoWaySplits] ...    
            = getStructureViolations(adjGraph);

% read RFC features from file

[s_skeletonEval.locallyConfidentSubgraphs,...
    s_skeletonEval.locallyWeakSubgraphs]...
    = subgraphEvaluation();

s_skeletonEval.abruptEnds = getAbruptEnds();



