function [A,b,senseArray,numEdges,numNodeConf,numRegions,nodeTypeStats]...
    = getILPConstraints(edgeListInds,edges2nodes,nodeEdgeIDs,junctionTypeListInds,...
        jEdges,dirEdges2regionsOnOff,setOfRegions,activeWSregionListInds_tr)

% version 3.1:
% 2014.01.07

% Inputs:
%   

% Outputs:
%   A - coefficient matrix for constraint equations. each row is a
%   constraint expression
%   columns of A: {eLIDs, eLIDs_complementary,eLIDs_off}{nodeConf}{regionsOn,regionsOff}
%       nodeConf = [{off,act1, act1',...}...]
%   b - RHS of constraint equations
%   senseArray - char array of =,<,> for each equation
%   numEdges - 
%   numNodeConf - number of active node configurations
%   numRegions - first region is the image border.
%   nodeTypeStats - 


%% Initialize

numEdges = numel(edgeListInds);
numRegions = size(setOfRegions,1) + 1; % + 1 for border

% node stats
% type 1 is J2 - junction with just 2 edges
[~, numJtypes] = size(jEdges);
% number of nodes for each type of junction
% col1: n.o nodes of this type
% col2: n.o tot junction activation variables \sum (nc2 * 2 + 1); +1 for
% offState
nodeTypeStats = zeros(numJtypes,3);
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
        nodeTypeStats(i,3) = numNodeConfVar_j * numJ_i; % num node confs for typeJ
    end
end

numNodeConf = sum(nodeTypeStats(:,3)); % tot. n.o. activation conf for all node types

numColsA = numEdges * 3 + numNodeConf + numRegions * 2;
% {edgeDir1}{edgeDir2}{edgeOff}{nodeConfigurations}{regionActive}{regionInactive}

% constraints types
% 1. Closedness with edge directionality (withCD)
% 2. Edge-region co-activation (withER)

withCD = 1;
withER = 1;
regionVarBin = 1;
enforceActiveRegions = 0;
enforceInactiveRegions = 0;

if(withCD)
    numCDeqns = numNodeConf;  % 
    disp('directionality & closedness constraint for node activation: ON ')
else
    numCDeqns = 0;
    disp('directionality & closedness constraint for node activation: OFF ')
end

if(withER)
    numEReqns = numEdges*2 + 1;
    disp('region and edge co-activation constraint: ON')
else
    numEReqns = 0;
    disp('region and edge co-activation constraint: OFF')
end

if(regionVarBin)
    numRVB = numRegions;
    disp('binary constraint for region variables: ON')
else
    numRVB = 0;
    disp('binary constraint for region variables: OFF')
end

if(enforceActiveRegions)
        numEnforceRegionsEqns = 1;
        disp('Enforce active regions from training labels: ON')
else
    numEnforceRegionsEqns = 0;
    disp('Enforce active regions from training labels: OFF')
end

if(enforceActiveRegions)
        numEnforceRegionInactEqns = 1;
        disp('Enforce inactive regions from training labels: ON')
else
    numEnforceRegionInactEqns = 0;
    disp('Enforce inactive regions from training labels: OFF')
end

totNumConstraints = numEdges + numNodeConf + numCDeqns + numEReqns ...
                    + numEnforceRegionsEqns + numEnforceRegionInactEqns...
                    + numRVB;

% Initialize outputs
A = sparse(totNumConstraints,numColsA);
b = zeros(totNumConstraints,1);
senseArray(1:totNumConstraints) = '=';

%% Equality constraint - edge activation
% Each edge in the graph correspond to 2 edge activation variables
% corresponding to the two directional edges between the pair of nodes
% connected by the original edge. Both components can't be active at the
% same time. The first direction corresponds to N1->N2 where N1 and N2 are
% given by the first and second cols of edges2nodes respectively.
% directional_edge_inds of of eLID = eLID, eLID+numEdges.
col_1 = 0;
for rowStop=1:numEdges
    col_1 = col_1 + 1;
    col_2 = col_1 + numEdges;
    col_3 = col_2 + numEdges;
    A(rowStop,col_1) = 1;
    A(rowStop,col_2) = 1;
    A(rowStop,col_3) = 1;
    b(rowStop) = 1;
    senseArray(rowStop) = '=';
end

%% Equality constraint node activation: 
% only one configuration per each node can be active the most.
% first nodeConf is the inactivation of the node
colStop = numEdges*3;
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

nodeConfStop = numEdges * 3; % points to the last filled node conf col ind in x
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
            nodeEdgeIDs_j(1) = [];  % first element is nodeInd
            nodeEdgeIDs_j = nodeEdgeIDs_j(nodeEdgeIDs_j>0);
            [~,nodeEdgeListInds_j] = intersect(edgeListInds',nodeEdgeIDs_j);
            % edgeListInds are the same as the col id for each edge
            
            % polarity (in/out) at each node
            edgePolarities_j = getEdgePolarity(edgeListInds,edges2nodes,nodeListInd_j,...
                    nodeEdgeListInds_j);
            % +1 for out. -1 for in
            
            % off state: all edges and all nodeConfs turned off (0)
            rowStop = rowStop + 1;
            b(rowStop) = 2;
            % senseArray(rowStop) = '='; default value
            nodeEdgeListInds_j_complementary = nodeEdgeListInds_j + numEdges;
            nodeEdgeListInds_j_dir = [nodeEdgeListInds_j; nodeEdgeListInds_j_complementary];
            
            nodeConfStop = nodeConfStop + 1;
            A(rowStop,nodeConfStop) = 2;
            A(rowStop,nodeEdgeListInds_j_dir) = 1;

            for k=1:numEdgeComb
                % for each edge combination
                edgePosInd = edgeCombinationList(k,:); % (1,2), (1,3) ..
                edgeLIDsInComb_k = nodeEdgeListInds_j(edgePosInd);
                polarities_k = edgePolarities_j(edgePosInd);
                % depending on the polarities make the two possible active
                % edge configurations where one edge comes in and the other
                % goes out
                edgeLID_1 = edgeLIDsInComb_k(1);
                edgeLID_2 = edgeLIDsInComb_k(2);
                complementaryEdgeLID_1 = numEdges + edgeLID_1;
                complementaryEdgeLID_2 = numEdges + edgeLID_2;
                sumPolarities_k = sum(polarities_k); 
                % sumPolarities_k gives an idea whether the edges are in
                % the same direction or opposing directions wrt to the
                % original direction assignment in edges2nodes.
                
                if(sumPolarities_k==2 || sumPolarities_k==-2)
                    % the 2 edges are in the same direction (in-in or out-out)
                    % assign opposing polarizations to the edges

                    % constraint for edge_1 of the pair: 1 0
                    % edge1p = 1
                    nodeConfStop = nodeConfStop + 1;

                    rowStop = rowStop + 1;
                    
                    A(rowStop,complementaryEdgeLID_1) = 1;
                    A(rowStop,edgeLID_2) = 1;
                    
                    A(rowStop,nodeConfStop) = -2;
                    b(rowStop) = -0.1;
                    senseArray(rowStop) = '>';
                    
                    rowStop = rowStop + 1;
                    A(rowStop,complementaryEdgeLID_1) = 1;
                    A(rowStop,edgeLID_2) = 1;
                    A(rowStop,nodeConfStop) = -2;
                    b(rowStop) = 1.1;
                    senseArray(rowStop) = '<';

                    % constraint for edge_2 of the pair 1 0
                    rowStop = rowStop + 1;
                    
                    A(rowStop,edgeLID_1) = 1;
                    A(rowStop,complementaryEdgeLID_2) = 1;
                    
                    nodeConfStop = nodeConfStop + 1;
                    A(rowStop,nodeConfStop) = -2;  % edge2p = 0
                    b(rowStop) = -0.1;
                    senseArray(rowStop) = '>';
                    
                    rowStop = rowStop + 1;
                    A(rowStop,edgeLID_1) = 1;
                    A(rowStop,complementaryEdgeLID_2) = 1;
                    A(rowStop,nodeConfStop) = -2;
                    b(rowStop) = 1.1;
                    senseArray(rowStop) = '<';

                elseif(sumPolarities_k==0)
                    % don't modify relative polarities

                    % 1. normal pair
                    rowStop = rowStop + 1;
                    
                    A(rowStop,edgeLID_1) = 1;
                    A(rowStop,edgeLID_2) = 1;
                    
                    nodeConfStop = nodeConfStop + 1;
                    A(rowStop,nodeConfStop) = -2;
                    b(rowStop) = -0.1;
                    senseArray(rowStop) = '>';
                    
                    rowStop = rowStop + 1;
                    A(rowStop,edgeLID_1) = 1;
                    A(rowStop,edgeLID_2) = 1;
                    A(rowStop,nodeConfStop) = -2;
                    b(rowStop) = 1.1;
                    senseArray(rowStop) = '<';
                    
                    % 2. complementary pair
                    rowStop = rowStop + 1;
                    
                    A(rowStop,complementaryEdgeLID_1) = 1;
                    A(rowStop,complementaryEdgeLID_2) = 1;
                    
                    nodeConfStop = nodeConfStop + 1;
                    A(rowStop,nodeConfStop) = -2;
                    b(rowStop) = -0.1; 
                    senseArray(rowStop) = '>';
                    
                    rowStop = rowStop + 1;
                    A(rowStop,complementaryEdgeLID_1) = 1;
                    A(rowStop,complementaryEdgeLID_2) = 1;
                    A(rowStop,nodeConfStop) = -2;
                    b(rowStop) = 1.1;
                    senseArray(rowStop) = '<';
                    
                else
                    disp('ERROR1:getILPConstraints.m problem with polarity calculation')
                    % disp('jType = %d   jInd = %d',jType,j)
                end
                
                % conf k1 (k)
                
                % conf k2 (k+1)
                
                
            end
           
       end
        
   end
end

%% Edge and region co-activation

if(withER)
    dirEdges2regionsOnOff = dirEdges2regionsOnOff + 1; 
    % regionID = 1 is for the image border, which should be off.
    r_offset = numEdges*3 + numNodeConf;
    % TODO: this for loop can be modified by considering the pair of
    % complementary edges inside the same loop. this can also be used to
    % check the accuracy of dirEdges2regionsOnOff
    for i=1:numEdges*2
        % each edge gives rise to 2 directional edge activation variables
        % find the two regions for the edge and if they're on or off
        rowStop = rowStop + 1;
        A(rowStop,i) = -1;
        
        % complementary edgeLID 
        if(i<=numEdges)
            edgeLID_comp = i + numEdges;
        else
            edgeLID_comp = i - numEdges;
        end
        A(rowStop,edgeLID_comp) = 1;
        
        rID_on = dirEdges2regionsOnOff(i,1);
        rID_off = dirEdges2regionsOnOff(i,2); 
        
        A(rowStop,(rID_on + r_offset)) = 1;
        A(rowStop,(rID_off + r_offset)) = -1;

        b(rowStop) = 0;
        % senseArray(rowStop) = '='; % default value
    end 
    % region 1 (image border) is off
    rowStop = rowStop + 1;
    rID_off = 1 + r_offset;
    A(rowStop,rID_off) = 1;
    b(rowStop) = 0;
    % senseArray(rowStop) = '='; % default value
end

%% Region variables should be binary
% rOn + rOff = 1;
if(regionVarBin)
    % each region var should be <=1
    colStart = numEdges*3 + numNodeConf + 1;
    colStop = numColsA - numRegions;
    
    for i=colStart:colStop
        rowStop = rowStop + 1;
        regionOnInd = i;
        regionOffInd = i+numRegions;
        A(rowStop,regionOnInd) = 1;
        A(rowStop,regionOffInd) = 1;
        b(rowStop) = 1;
        senseArray(rowStop) = '=';
    end
    
end
%% Enforce Active regions
if(enforceActiveRegions)
    
    numActiveRegions_tr = numel(activeWSregionListInds_tr);
    r_offset = numEdges*3 + numNodeConf;
    offset_activeRegionWSInds = r_offset + activeWSregionListInds_tr;
    rowStop = rowStop + 1;
    A(rowStop,offset_activeRegionWSInds) = 1;
    b(rowStop) = numActiveRegions_tr - 10;
    senseArray(rowStop) = '>';
end

if(enforceInactiveRegions)
   regionSeq = 1:numRegions;
   inactiveWsRegionListInds_tr = setdiff(regionSeq,activeWSregionListInds_tr);
   r_offset = numEdges*3 + numNodeConf;
   offset_inactiveRegionWsInds = r_offset + inactiveWsRegionListInds_tr;
   
   rowStop = rowStop + 1;
   A(rowStop,offset_inactiveRegionWsInds) = 1;
   b(rowStop) = 10;
   senseArray(rowStop) = '<';
end

