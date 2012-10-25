%% function IOut = sampleFromMRF_simple(imgDim,bb,Hmat,Vmat,prior,Dictionary)

% generates a sample image based on the learned MRF. Generates the image
% sequentially from the top left corner to the bottom right corner.

% Inputs:
% imgDim - vector containing the size of the output image to be sampled [rows columns]
% bb - block size
% Hmat - horizontal association matrix learned (row,column -> left to right)
% Vmat - vertical association matrix learned (row,column -> top to bottom)
% prior - prior for label (word) usage
% Dictionary

% Output:
% Iout - sampled image

function IOut = sampleFromMRF_simple(imgDim,bb,Hmat,Vmat,Dictionary) 

                                        % Initialize random field
rows = imgDim(1) - bb + 1;              % number of overlapping image patches per column
cols = imgDim(2) - bb + 1;              % number of overlapping image patches per row

totPatches = rows*cols;                 % total number of patches to be generated
coefMat = sparse(size(Dictionary,2),totPatches);

%% picking appropriate words for each overlapping patch of the image
% Pick a random word for the first image patch
wordInd = ceil(rand(1) * size(Dictionary,2)) ; % use prior ? TODO
coefMat(wordInd,1) = 1;                 % assign the word to the first patch
for i = 2:totPatches
                                        % loop until all the patches are processed
    if(mod(i,cols) ~= 1) % i is not in the first colum
        wordInFirstCol = 0;             % flag
                                        % sample a word from Hmat given i-1 th word
        wordToLeftInd = find(coefMat(:,i-1));                                
                                        
        [wordIndH,probH] = sampleWordToRight(wordToLeftInd,Hmat);    
    else
        wordInFirstCol = 1;
    end
    if(i>cols)                          % i is not in the first row
        wordInFirstRow = 0;
                                        % sample a word from Vmat given the one above
        wordToTopInd = find(coefMat(:,i-cols));                                
        [wordIndV,probV] = sampleWordToBottom(wordToTopInd,Vmat);  
    else
        wordInFirstRow = 1;
    end
        
    % assigning a word to the patch i based on the above information
    if(wordInFirstCol)
        % we only have to consider the vertical association
        wordInd = wordIndV;
    elseif (wordInFirstRow)
        % we only have to consider the horizontal associations
        wordInd = wordIndH;
    else
        % we have to consider both associations
        % pick the stronger candidate based on the probability
        if (probH > probV)
            wordInd = wordIndH;
        elseif (probH == probV)
            if(rand(1)<0.5)
                wordInd = wordIndH;
            else
                wordInd = wordIndV;
            end
        else
            wordInd = wordIndV;
        end
        
    end
    
    coefMat(wordInd,i) = 1;             % assigning the word to the patch i
    
end

%% Combine the overlapping patches to get the sampled whole image
IOut = getImageFromCoeff(Dictionary,coefMat,imgDim);
figure(3)
imagesc(IOut);

