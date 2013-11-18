function [A,b,senseArray,numEdges,numNodeConf,numRegions,nodeTypeStats]...
    = getILPConstraints(edgeListInds,edges2nodes,nodeEdgeIDs,junctionTypeListInds,...
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
% col2: n.o tot junction activation variables \sum (nc2 * 2 + 1); +1 for
% offState
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
        numNodeConfVar_j = numEdgeCombinations *2 + 1; % n.o conf per node
        nodeTypeStats(i,2) = numNodeConfVar_j;
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
% first nodeConf is the inactivation of the node
colStop = numEdges*2;
rowStop = 0;
for jType = 1:numJtypes
    % for each junction type
    numNodes_j = nodeTypeStats(jType,1);
    numCoef = nodeTypeStats(jType,2); % n.o nodeState var per node, incl offState
    rowStart = rowStop + 1; 
    rowStop = rowStop + numNodes_j;
    for row=rowStart:rowStop
        colStart = colStop + 1;
        colStop = colStop + numCoef;
        A(row,colStart:colStop) = 1;
        b(row) = 1;
        % senseArray already contains the correct symbol '='
    end    
end

%% Equallity constraint: edge directionality at active nodes 
% each node configuration is an active configuration where exactly 2 edges
% have to be active.
% the polarity of the edge is set to 0 if the direction implied in
% edges2nodes is compliant with
nodeConfStop = numEdges * 2; % points to the last filled node conf col ind in x
if(withCD)

   for jType=1:numJtypes
       numJ_i = nodeTypeStats(jType,1);
       numEdges_i = jType +1;
       numEdgeComb = nchoosek(numEdges_i,2);
       numNodeConf_j = numEdgeComb * 2;
       edgeSeq = 1:numEdges_i;
       edgeCombinationList = nchoosek(edgeSeq,2);
       for j=1:numJ_i
            % what edges are connected to this node
            nodeListInd_j = junctionTypeListInds(j,jType);
            nodeEdgeIDs_j = nodeEdgeIDs(nodeListInd_j,:);
            nodeEdgeIDs_j = nodeEdgeIDs_j(nodeEdgeIDs_j>0);
            [~,nodeEdgeListInds_j] = intersect(edgeListInds,nodeEdgeIDs_j);
            % edgeListInds are the same as the col id for each edge
            
            % polarity (in/out) at each node
            edgePolarities_j = getEdgePolarity(edgeListInds,edges2nodes,nodeListInd_j);
            

            % off state: all edges and all nodeConfs turned off (0)
            rowStop = rowStop + 1;
            b(rowStop) = 2;
            % senseArray(rowStop) = '='; default value
            nodeConfStop = nodeConfStop + 1;
            A(rowStop,nodeConfStop) = 2;
            A(rowStop,nodeEdgeListInds_j) = 1;

            for k=1:numEdgeComb
                % for each edge combination
                edgePosInd = edgeCombinationList(k,:); % (1,2), (1,3) ..
                edgeLIDsInComb_k = nodeEdgeListInds_j(edgePosInd);
                polarities_k = edgePolarities_j(edgeLIDsInComb_k);
                % depending on the polarities make the two possible active
                % edge configurations where one edge comes in and the other
                % goes out
                edgeLID_1 = edgeLIDsInComb_k(1);
                edgeLID_2 = edgeLIDsInComb_k(2);
                sumPolarities_k = sum(polarities_k);
                rowStop = rowStop + 1;
                nodeConfStop = nodeConfStop + 1;
                if(sumPolarities_k==2 || sumPolarities_k==0)
                    % assign opposing polarizations to the edges
                elseif(sumPolarities_k==1)
                    % assign the same polarization to the edges
                else
                    disp('ERROR1:getILPConstraints.m problem with polarity calculation')
                    disp('jType = %d   jInd = %d',jType,j)
                end
                
                % conf k1 (k)
                
                % conf k2 (k+1)
                
                
            end
           
       end
        
   end
end

%% Edge and region co-activation


