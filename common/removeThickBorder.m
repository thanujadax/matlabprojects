function outImg = removeThickBorder(inImg,margin)

[inSizeR, inSizeC, dim] = size(inImg);

outSizeR = inSizeR - margin*2;
outSizeC = inSizeC - margin*2;

outImg = zeros(outSizeR,outSizeC,dim);

startR = margin + 1;
stopR = inSizeR - margin;

startC = margin + 1;
stopC = inSizeC - margin;

% m=0;
% n=0;
% 
% for k=1:dim
%     for i=startR:stopR
%         m = m +1;
%         for j=startC:stopC
%             n = n + 1;
%             outImg(m,n,dim) = inImg(i,j,dim); 
%         end
%     end
% end


for k=1:dim
    inImg_k = zeros(inSizeR,inSizeC);
    outImg_k = zeros(outSizeR,outSizeC);
    inImg_k = inImg(:,:,k);
    outImg_k = inImg_k(startR:stopR,startC:stopC);
    outImg(:,:,k) = outImg_k; 
end