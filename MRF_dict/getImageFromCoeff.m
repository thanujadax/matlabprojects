function IOut = getImageFromCoeff(Dictionary,coefMat,imgSize)

NN1 = imgSize(1);
NN2 = imgSize(2);

blocks = Dictionary*coefMat;        % blocks have image patches as its columns
numBlocks = (NN1-bb+1)*(NN2-bb+1);
idx = 1:numBlocks;
count = 1;
Weight = zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(imgSize-bb+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
    block = reshape(blocks(:,count),[bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
    count = count+1;
end

IOut = IMout./Weight;