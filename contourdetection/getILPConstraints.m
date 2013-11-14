function [A,b,senseArray,numEdges,numNodeConf,numRegions]...
    = getILPConstraints(edgeListInds,edges2nodes,nodeEdges,junctionTypeListInds,...
        jEdges,twoRegionEdges,edges2regions,setOfRegionsMat)

% version 3:
% 2013 11 12

% Inputs:
%   

% Outputs:
%   A - coefficient matrix for constraint equations. each row is a
%   constraint expression
%   b - RHS of constraint equations
%   senseArray - char array of =,<,> for each equation
%   numEdges - 
%   numNodeConf - number of active node configurations
%   numRegions - 

%% Initialize

numEdges = numel(edgeListInds);
numRegions = size(setOfRegionsMat,1);

% node stats
% type 1 is J2 - junction with just 2 edges
[~, numJtypes] = size(jEdges);
% number of nodes for each type of junction
% col1: n.o nodes of this type
% col2: n.o tot junction activation variables \sum nc2 * 2
nodeTypeStats = zeros(numJtypes,2);
% 
for i=1:numJtypes
    jEdges_i = jEdges{i};
    if(jEdges_i==0)
        % no junctions of this type
        nodeTypeStats(i,1) = 0; % n.o nodes
        nodeTypeStats(i,2) = 0; % n.o tot junction activation variables
    else
        [numJ_i,numEdges_i] = size(jEdges_i);
        nodeTypeStats(i,1) = numJ_i;
        
        numEdgeCombinations = nchoosek(numEdges_i,2);
        numActivationVar = numEdgeCombinations * 2 * numJ_i;
        nodeTypeStats(i,2) = numActivationVar;
    end
end

numNodeConf = sum(nodeTypeStats(:,2));

numColsA = numEdges * 2 + numNodeConf + numRegions;




%% Equallity constraint: node activation with edge directionality

