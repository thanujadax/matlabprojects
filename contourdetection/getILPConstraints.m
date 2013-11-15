function [A,b,senseArray,numEdges,numNodeConf,numRegions,nodeTypeStats]...
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
        nodeTypeStats(i,1) = numJ_i; % n.o nodes
        
        numEdgeCombinations = nchoosek(numEdges_i,2); 
        numActivationVar = numEdgeCombinations * 2 * numJ_i; % n.o active conf in tot_i
        nodeTypeStats(i,2) = numActivationVar;
    end
end

numNodeConf = sum(nodeTypeStats(:,2)); % tot. n.o. activation conf for all node types

numColsA = numEdges * 2 + numNodeConf + numRegions;
% {edgeActivation}{edgePolarity}{nodeConfigurations}{regionActivation}

% constraints types
% 1. Closedness with edge directionality (withCD)
% 2. Edge-region co-activation (withER)

withCD = 1;
withER = 1;

if(withCD)
    numCDeqns = numNodeConf;  % 
    disp('directionality & closedness constraint for node activation: ON ')
else
    numCDeqns = 0;
    disp('directionality & closedness constraint for node activation: OFF ')
end

if(withER)
    numEReqns = numRegions;
    disp('region and edge co-activation constraint: ON')
else
    numEReqns = 0;
    disp('region and edge co-activation constraint: OFF')
end


totNumConstraints = numCDeqns + numEReqns;

% Initialize outputs
A = zeros(totNumConstraints,numColsA);
b = zeros(totNumConstraints,1);
senseArray(1:totNumConstraints) = '=';

%% Equality constraint node activation: 
% only one configuration per each node can be active the most.
colStop = 0;
rowStop = 0;
for jType = 1:numJtypes
    % for each junction type
    numNodes_j = nodeTypeStats(jType,1);
    numEdgePJ = jType + 1;      % number of edges per junction
    numCoef = nchoosek(numEdgePJ,2)*2 + 1; % num edge pair combinations + inactivation  
    rowStart = rowStop + 1; 
    rowStop = rowStart - 1 + numNodes_j;
    for row=rowStart:rowStop
        colStart = colStop + 1;
        colStop = colStart - 1 + numCoef;
        A(row,colStart:colStop) = 1;
        b(row) = 1;
        % senseArray already contains the correct symbol '='
    end    
end

%% Equallity constraint: edge directionality at active nodes 

if(withCD)
   % directionality and closedness constraint for node activation
   for jType=1:numJtypes
       numJ_i = nodeTypeStats(jType,1);
       numEdges_i = jType +1;
       numEdgeComb = nchoosek(numEdges_i,2);
       for j=1:numJ_i
           
       end
        
   end
end

