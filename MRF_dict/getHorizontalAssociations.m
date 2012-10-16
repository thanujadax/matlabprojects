function AMat = getHorizontalAssociations(coefmat)

% coefmat - sparse matrix of coefficients. Rows correspond to the
% dictionary words and the columns correspond to the overlapping patches of
% the image.

% % fixed parameters
% numwords = 1;           % number of active coefficients per patch
% HorizontalAssociations = ones(size(Dictionary));
% 
% 
% [NN1,NN2] = size(image);
% 
% while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
%     slidingDis = slidingDis+1;
% end
% [blocks,idx] = my_im2col(image,[bb,bb],slidingDis);
% 
% % processing 30000 blocks at a time
% for jj = 1:30000:size(blocks,2)
%     jumpSize = min(jj+30000-1);
%     
%     if(reduceDC)
%         vecOfMeans = mean(blocks(:,jj:jumpSize));
%         blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
%     end
%     
%     Coefs = OMP([FixedDictionaryElement,Dictionary],blocks, numwords);
%     if(reduceDC)
%         blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
%     else
%         blocks(:,jj:jumpSize)= Dictionary*Coefs ;
%     end
%     
% end

AMat = ones(size(coefmat,1));     % initializing association matrix

ind1 = find(coefmat(:,1));       % index of the nonzero element of column 1 (index of atom for patch #1)
for i = 2:size(coefmat,2)
    ind2 = find(coefmat(:,i));  % index of the word for the next patch
    AMat(ind1,ind2) = AMat(ind1,ind2) + 1;
    ind = ind2;
end


