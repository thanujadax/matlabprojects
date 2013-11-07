function features = writeLabelsFile(labelVector)

% labels.txt


filename = 'labels.txt';
fileID = fopen(filename,'w');

numVar = numel(labelVector);

for i=1:numVar
    ftval = '%d \n';
    fprintf(fileID,ftval,labelVector(i));
end

fclose(fileID);
features = 0;