function featureMat = writeFeaturesFile_3DILP(unaryScoresMat,pathToSave)

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

[numGridCells,numFeatures] = size(unaryScoresMat);

numRows = numGridCells * numFeatures;

featureMat = zeros(numRows,numFeatures);

filename = 'features.txt';
% filename = fullfile(pathToSave,filename);
fileID = fopen(filename,'w');

rowStop = 0;
for i=1:numFeatures
    rowStart = rowStop + 1;
    rowStop = rowStop + numGridCells;
    featureMat(rowStart:rowStop,i) = unaryScoresMat(:,i);
end


%%  write to file
for i=1:numRows
   fprintf(fileID, '%4.6f ', featureMat(i,:));
   fprintf(fileID, '\n');
end
fclose(fileID);
