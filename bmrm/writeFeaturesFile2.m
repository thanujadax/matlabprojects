function featureMat = writeFeaturesFile2(f,jEdges,numEdges,numRegions,outputPath)

% features.txt
% 
% # contains the feature vectors for the variables, one per row
% #
% # x x x x x 0 0 0 0 0 0 0 0
% # x x x x x 0 0 0 0 0 0 0 0
% # x x x x x 0 0 0 0 0 0 0 0
% # 0 0 0 0 0 x x x x x 0 0 0 |
% # 0 0 0 0 0 x x x x x 0 0 0 | one category of variables
% # 0 0 0 0 0 x x x x x 0 0 0 |
% # 0 0 0 0 0 x x x x x 0 0 0 |
% # 0 0 0 0 0 0 0 0 0 0 x x x
% # 0 0 0 0 0 0 0 0 0 0 x x x
% # 0 0 0 0 0 0 0 0 0 0 x x x
% #
% # different training sets can just be concatenated

numFeatures = 6;
numRows = numel(f);

featureMat = zeros(numRows,numFeatures);

filename = 'features.txt';
filename = fullfile(outputPath,filename);
fileID = fopen(filename,'w');

numVar = numel(f);

% for i=1:numVar
%     ftval = '%3.5f \n';
%     fprintf(fileID,ftval,f(i));
% end

% fclose(fileID);
% features = 0;


%% edges
% edge on dir1 & dir2
featInd = 0;
for i=1:numEdges*2
    featInd = i;
    featVal = f(featInd);  
    featureMat(featInd,1) = featVal;  % (1) w_e_on     
end
% edge off
featInd = numEdges * 2;
for i=1:numEdges
    featInd = featInd + 1;
    featVal = f(featInd); 
    featureMat(featInd,2) = featVal;  % (2) w_e_off 
end

%% nodes

% create label vector
[~, numJtypes] = size(jEdges);
% type 1 is J2 - junction with just 2 edges
nodeTypeStats = zeros(numJtypes,2);
% each row corresponds to a junction type. row 1: type 1 (J2)
% column 1: n.o. junction nodes of this type
% column 2: n.o. number of node configurations (active and inactive)
totJunctionVar = zeros(numJtypes,1); % stores the number of coefficients for each type of J
totActiveJunctionConfs = zeros(numJtypes,1);

for i=1:numJtypes
    jEdges_i = jEdges{i};
    if(jEdges_i==0)
        % no junctions of this type
        nodeTypeStats(i,1) = 0;
        nodeTypeStats(i,2) = 0;
        totJunctionVar(i) = 0;
        totActiveJunctionConfs(i) = 0;
    else
        [numJ_i,numEdges_i] = size(jEdges_i);
        nodeTypeStats(i,1) = numJ_i;
        
        numEdgeCombinations = nchoosek(numEdges_i,2);
        nodeTypeStats(i,2) = numEdgeCombinations + 1; % +1 for the inactivation state
        
        totJunctionVar(i) = numJ_i.*(numEdgeCombinations + 1); % 1 for the inactive state
        totActiveJunctionConfs(i) = numJ_i.*numEdgeCombinations;
        % clear nodeAngleCost_i
    end
end

for i=1:numJtypes
    % for each junction type
    clear nodeAngleCost_i
    nodeEdgeIDs_i = jEdges{i};
    numNodes_i = size(nodeEdgeIDs_i,1);
    numCoeff_i = nodeTypeStats(i,2);
    numActiveStates = numCoeff_i - 1;
    
    if( sum(sum(~isnan(nodeEdgeIDs_i)))>0 && nodeEdgeIDs_i(1)~=0)
        for j=1:numNodes_i
            featInd = featInd + 1;
            featVal = f(featInd);
            featureMat(featInd,3) = featVal; % (3) w_n_off
            
            featStart = featInd + 1;
            featInd = featInd + numActiveStates;
            
            featureMat(featStart:featInd,4) = f(featStart:featInd); % (4) w_n_on
            
            % w_n_on complementary directions of edges (same feature vales)
            featStart = featInd + 1;
            featInd = featInd + numActiveStates;
            
            featureMat(featStart:featInd,4) = f(featStart:featInd); % (4) w_n_on
            
        end
    end
end


%% regions
% rOn
for i=1:numRegions
    featInd = featInd + 1;
    featVal = f(featInd);
    featureMat(featInd,5) = featVal;  % (5) w_r_on
end
% rOff
for i=1:numRegions
    featInd = featInd + 1;
    featVal = f(featInd);
    featureMat(featInd,6) = featVal;  % (6) w_r_off
end

%%  write to file
for i=1:numVar
   fprintf(fileID, '%4.6f ', featureMat(i,:));
   fprintf(fileID, '\n');
end
fclose(fileID);
