function writeImageMatrixToTiffStack(imageMatrix,outputFileName)
% remove outputfile if already exists
if exist(outputFileName,'file') == 2
    delete(outputFileName);
end
    
    
for K=1:length(imageMatrix(1, 1, :))
   imwrite(imageMatrix(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
end