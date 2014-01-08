function f = getILPObjectiveVectorParametric2(edgeUnary,nodeAngleCosts,...
            regionUnary,w_on_e,w_off_e,w_off_n,w_on_n,w_on_r,w_off_r,...
            nodeTypeStats)
        
% version 2 with directed edges - 2014.01.07
% {edges}{edges_complementary}{nodeConfigurations}{border,regionActivation}

% Inputs:
%   nodeTypeStats: each row corresponds to a junctionType (n.o edges connected)
%       col1: number of nodes of  this junctionType
%       col2: total number of junction variables = \sum (nc2 * 2 + 1); +1 for this type 
%   edgeUnary: unary cost for directed edge activation in the order given
%   in edgesListInds_directional


numEdges = numel(edgeUnary);
numEdges_directed = numEdges * 2;

[~, numJtypes] = size(nodeAngleCosts);
% type 1 is J2 - junction with just 2 edges

totJunctionVar = sum(nodeTypeStats(:,3));
numRegions = numel(regionUnary) + 1; % region 1 is image border
numElements = numEdges*3 + totJunctionVar + numRegions*2;
f = zeros(numElements,1);
% order of elements in f
%[{edges ... edges_comp ... edgeOff },{j01,j11,..,j02,j12,...},{rOn... rOff}]

%% edge variables
% for each edge, there are 2 possible active states. Only one of those can
% be maximally active per each edge (see constraints file). However,
% different costs should be allocated for the two variables since one
% direction (N1->N2) is more preferred by the data.
% e.g. for starters, we can assign the negative of the unary value to the
% less preferred direction

% j=1;
% for i=1:2:2*numEdges
%     f(i) = edgeUnary(j) * w_off_e;        % inactivation cost *w_off
%     f(i+1) = edgeUnary(j) * w_on_e;     % activation cost for the same edge *w_on
%     j = j+1;
% end

for i=1:numEdges
    f(i) = edgeUnary(i) * w_on_e;
    f(i+numEdges) = edgeUnary(i) * w_on_e;
    f(i+numEdges*2) = edgeUnary(i) * w_off_e;
%     
end

f_stop_ind = numEdges*3;

%%  junction variables
for i=1:numJtypes
    % for each junction type
    clear nodeAngleCost_i
    nodeAngleCost_i = nodeAngleCosts{i};
    notNaN = sum(sum(~isnan(nodeAngleCost_i)))>0;
    if(~isempty(nodeAngleCost_i) && notNaN)
        [numJ,numc] = size(nodeAngleCost_i);
        nodeAngleCost_bidir_i = zeros(numJ,(numc*2));
        nodeAngleCost_bidir_i = nodeAngleCost_bidir_i .* w_on_n;

        % inactivation cost - comes from learned param
        offJcost = ones(numJ,1).*w_off_n;

        % REMOVED: bbJunction cost assignment.

        nodeAngleCost_bidir_i = [offJcost nodeAngleCost_bidir_i];
        numCoeff_i = nodeTypeStats(i,3);
        f_start_ind = f_stop_ind + 1;
        f_stop_ind = f_start_ind + numCoeff_i - 1;
        % assign coefficients to vector f
        angleCostMat_i = nodeAngleCost_bidir_i';
        f(f_start_ind:f_stop_ind) = angleCostMat_i(1:numCoeff_i);
    end
end

%% cells
k=1;
% regionID = 1 is the image border. this is inactive (see constraints file)
f_stop_ind = f_stop_ind + 1;
f(f_stop_ind) = 0; % for the border

f_stop_ind = f_stop_ind + 1;
for i=(f_stop_ind):(f_stop_ind+numRegions-2)
    f(i) = regionUnary(k) * w_on_r; % activation
    f(i+numRegions) = regionUnary(k) * w_off_r; % inactivation
    k = k + 1;
end