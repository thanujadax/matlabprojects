% function IOut = OMPDenoisedImage(Image,Dictionary0,DictionaryN,bb,...
%     maxBlocksToConsidersigma,C,slidingDis,waitBarOn,reduceDC)
% 
% Takes in the image to be denoised, the first dictionary and the second
% dictionary and produces a denoised image  using the words of the second
% dictionary, based on the coefficients learned via  the first dictionary
%
% Inputs:
% Image - image to be denoised
% Dictionary0 - denoised dictionary
% DictionaryN - initial dictionary
% bb - block size
% maxBlocksToConsider - to adjust the sliding distance (for defining blocks)
% sigma
% C
% slidingDis
% waitBarOn
% numWords: if zero, go for the error goal. If nonzero, this is the number
% of atoms allowed per represented patch.
%
% Thanuja 


function [IOut, allcoeffs,vecOfMeans] = SparseCodeImageNN(Image,Dictionary,bb,...
    maxBlocksToConsider,sigma,C,slidingDis,waitBarOn,reduceDC,numWords,numIter)

vecOfMeans = 0;

[NN1,NN2] = size(Image);
errT = sigma*C;      % Error threshold for OMP (OMPErr)

%blocks = im2col(Image,[NN1,NN2],[bb,bb],'sliding');
while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end
[blocks,idx] = my_im2col(Image,[bb,bb],slidingDis);

if (waitBarOn)
    % newCounterForWaitBar = (param.numIteration+1)*size(blocks,2);
end

% Coefs = sparse(param.K,size(blocks,2));
%Coefs = sparse(size(Dictionary,2),stepWidth);
allcoeffs = sparse(0,0);
% sparsecoeff = zeros(size(Dictionary,2),size(blocks,2));
% sparsecoeff = 0;
% go with jumps of stepWidth

stepWidth=100;

% calculating coefficients with the non-negativity constraint
for jj = 1:stepWidth:size(blocks,2)
    jumpSize = min(jj+stepWidth-1,size(blocks,2));
    if (reduceDC)
        vecOfMeans = mean(blocks(:,jj:jumpSize));
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
    end

    display('processing patch ');
    display(jj);
    display(' to ');
    display(jumpSize);
    
    if(numWords==1)
        Coefs = least_squares_match(blocks(:,jj:jumpSize), Dictionary);
        Coeff = full(Coefs);
        %Coefs = OMP(Dictionary,blocks(:,jj:jumpSize),numWords);
    else
        Coeff = OMPerr(Dictionary,blocks(:,jj:jumpSize),errT);
    end
    
    if (reduceDC)
        blocks(:,jj:jumpSize)= Dictionary*Coeff + ones(size(blocks,1),1) * vecOfMeans;
    else
        blocks(:,jj:jumpSize)= Dictionary*Coeff ;
    end
    
    % kk = size(Coeff,2)+jj-1;
    % sparsecoeff(:,jj:kk) = Coefs;
    allcoeffs = [allcoeffs Coefs];

end


% for jj = 1:stepWidth:size(blocks,2)
%     if (waitBarOn)
%         % waitbar(((param.numIteration*size(blocks,2))+jj)/newCounterForWaitBar);
%     end
%     jumpSize = min(jj+stepWidth-1,size(blocks,2));
%     if (reduceDC)
%         vecOfMeans = mean(blocks(:,jj:jumpSize));
%         blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
%     end
%     
%     % Coefs = mexOMPerrIterative(blocks(:,jj:jumpSize),Dictionary,errT);
%     %Coefs = OMPerrOne(Dictionary,blocks(:,jj:jumpSize),errT);
%     if(numWords>0)
%         Coefs = OMP(Dictionary,blocks(:,jj:jumpSize),numWords);
%     else
%         Coefs = OMPerr(Dictionary,blocks(:,jj:jumpSize),errT);
%     end
%     
%     if (reduceDC)
%         blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
%     else
%         blocks(:,jj:jumpSize)= Dictionary*Coefs ;
%     end
%     
%     kk = size(Coefs,2)+jj-1;
% %     if(kk>size(sparsecoeff,2))
% %         kk = size(sparsecoeff,2);
% %     end
%     sparsecoeff(:,jj:kk) = Coefs;
% end

count = 1;
Weight = zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(Image)-bb+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
    block = reshape(blocks(:,count),[bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
    count = count+1;
end;

if (waitBarOn)
%    close(h);
end
%IOut = (Image+0.034*sigma*IMout)./(1+0.034*sigma*Weight);
IOut = IMout./Weight;