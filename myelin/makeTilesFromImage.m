function makeTilesFromImage(inputImageFileName,params,outputSavePath)

% params.xBorder = 25; % px
% params.yBorder = 25; % px
% params.xTiles = 40;
% params.yTiles = 40;

inputImage = imread(inputImageFileName);
[sizeR,sizeC] = size(inputImage);

tileRoISizeC = floor((sizeC-2*params.xBorder)/params.xTiles);
tileRoISizeR = floor((sizeR-2*params.yBorder)/params.yTiles);

for i = 1:params.xTiles
    tileCstart = (i-1)*tileRoISizeC + 1;
    tileCstop = tileCstart + tileRoISizeC +2*params.xBorder -1;
    for j = 1:params.yTiles
        tileRstart = (j-1)*tileRoISizeR + 1;
        tileRstop = tileRstart + tileRoISizeR + 2*params.yBorder - 1;
        
        tile = inputImage(tileRstart:tileRstop,tileCstart:tileCstop);
        tokenizedFileName = strsplit(inputImageFileName,filesep);
        fileName = tokenizedFileName{end};
        fileName = strtok(fileName,'.');
        fileName = sprintf('%s_tinyTile_row%02d_col%02d.png',fileName,j,i);
        saveFileName = fullfile(outputSavePath,fileName);
        imwrite(tile,saveFileName,'png');
    end
end