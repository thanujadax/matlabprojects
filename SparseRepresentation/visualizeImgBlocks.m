function imgBlocksAsImg = visualizeImgBlocks(inputImage,bb,slidingDist)

% extract overlapping image patches as column vectors
[imgblocks,idx] = my_im2col(inputImage,[bb,bb],slidingDist);

[dimX dimY] = size(inputImage);

numBlksPerRow = dimX - bb + 1;
numBlksPerCol = dimY - bb + 1;

patchIndex = 1;
for i = 1:numBlksPerRow
    startingPoint_i = (i-1)*bb + 1;
    for j = 1:numBlksPerCol
        startingPoint_j = (j-1)*bb + 1;
        imgBlocksAsImg(startingPoint_i:(startingPoint_i+bb-1),...
            startingPoint_j:(startingPoint_j+bb-1)) = ...
                reshape(Dictionary(:,labelVector(patchIndex)),bb,bb);
        patchIndex = patchIndex + 1;
        
    end
end



