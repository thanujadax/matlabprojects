function features = writeLabelsFile2(labelVector,pathToSave)

% labels.txt


filename = 'labels.txt';
% filename = fullfile(pathToSave,filename);
fileID = fopen(filename,'wt');
fprintf(fileID,'%d\n',labelVector);
fclose(fileID);
features = 0;

