% create one tiff stack out of many tiff stacks
function combineTifPagesToSingleTif()

inputDir = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/em_2013january/raw/';
outputFileName = '/home/thanuja/projects/data/drosophilaLarva_ssTEM/rawStack.tif';

inputContent = dir(fullfile(inputDir,'*.tif'));

   fileName = inputContent(1).name;
   fileName = fullfile(inputDir,fileName);
   image = readTiffStackToArray(fileName);
   
[sizeR,sizeC] = size(image);

syntheticStack = zeros(sizeR,sizeC,length(inputContent));

for i=1:length(inputContent)
    
   fileName = inputContent(i).name;
   fileName = fullfile(inputDir,fileName);
   image = readTiffStackToArray(fileName);
    
   syntheticStack(:,:,i) = image(:,:);
end

syntheticStack = syntheticStack./255;

for K=1:size(syntheticStack,3)
    imwrite(syntheticStack(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
end