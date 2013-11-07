function features = writeFeaturesFile(f)

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

filename = 'features.txt';
fileID = fopen(filename,'w');

numVar = numel(f);

for i=1:numVar
    ftval = '%3.5f \n';
    fprintf(fileID,ftval,f(i));
end

fclose(fileID);
features = 0;
