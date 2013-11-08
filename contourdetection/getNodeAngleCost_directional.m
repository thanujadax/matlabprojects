function nodeAngleCost = getNodeAngleCost_directional(theta,alpha,edgePriors_j,...
                        w_on_n)
                    
% calculate the cost for each active configuration of a node given the edge
% angles according to the OFR (theta) and the edge position of the graph
% relative the current node (alpha)

% the number of columns of each theta and alpha corresponds to the number
% of edges connected to a node. this gives rise to a combination of the
% pairs of edges we can pick - nchoosek(n,2). in the output: nodeAngleCost,
% we have a number of columns equal to the number of such combinations

% edgePriors_j contains the edge priors for all the nodes in the junction
% type in concern, in the order considered in theta and alpha


% 20131108: assigned the same weight w_on_n to both negative and positive costs


% Parameters
K = 1;    % scalar factor for node cost 

if(theta==0)
    nodeAngleCost = nan;
elseif(~isempty(theta))
    [numNodes,numEdgesPerNode] = size(theta);
    if(numEdgesPerNode>0)
        numCombinations = nchoosek(numEdgesPerNode,2);
        edgeIDvect = 1:numEdgesPerNode;
        combinations = nchoosek(edgeIDvect,2);

        % get the outwardness score for all the edge for all the nodes
        outwardnessScores = getOutwardness(theta,alpha);

        % calculate the cost for each pair of edge combinations - nodeAngleCost
        % = multiplication of the outwardness scores of the two edges
        nodeAngleCost = zeros(numNodes,numCombinations);
        for i=1:numNodes
           for j=1:numCombinations
                edge1LInd = combinations(j,1);
                edge2LInd = combinations(j,2);
                edgePriorFactor = edgePriors_j(i,edge1LInd) * edgePriors_j(i,edge2LInd)*K;
                nodeAngleCost(i,j) = outwardnessScores(i,edge1LInd) *...
                                    outwardnessScores(i,edge2LInd) * edgePriorFactor;

           end
        end
        ind_neg = (nodeAngleCost<0);
        if(~isempty(ind_neg))
            negCosts = nodeAngleCost(ind_neg);
            % negCosts = negCosts .* cNeg;    % further scaling the negative values (only)
            negCosts = negCosts .* w_on_n;
            nodeAngleCost(ind_neg) = negCosts;
        end

        ind_pos = (nodeAngleCost>0);
        if(~isempty(ind_pos))
            posCosts = nodeAngleCost(ind_pos);
%             posCosts = posCosts .* cPos;    % further scaling of the positive values (only)
            posCosts = posCosts .* w_on_n;
            nodeAngleCost(ind_pos) = posCosts;
        end
    else
        nodeAngleCost = 0;
    end
else
    nodeAngleCost = 0;
end