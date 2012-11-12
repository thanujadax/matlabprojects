function Coefs = XwordPick(Data,Dictionary,bb,invertImages)

% selects the best matching Dictionary word per each image patch in Data
% based on correlation. It's essential that the structures we look for have
% higher values in the inputs. i.e. EM images have to be inverted.

% inputs:
% Data - input image blocks as a set of column vectors
% Dictionary - 
% invertImages - if 1, the images should be inverted

if(invertImages)
    Data = invertImage(Data);
    Dictionary = invertImage(Dictionary);
end
    
Coefs = sparse(size(Dictionary,2),size(Data,2));

for i = 1:size(Data,2)
    % for each input image patch i
    % calculate correlations
    dataPatch = reshape(Data(:,i),bb,bb);
    bestCorr = -1;
    bestWord = -1;

    for j = 1:size(Dictionary,2)
       % for each Dictionary word j 
       % to be correlated with the i-th image patch in data
       dictionaryPatch = reshape(Dictionary(:,j),bb,bb); 
       corr = xcorr2(dataPatch,dictionaryPatch);
       if(max(max(corr))>bestCorr)
           bestWord = j;
           bestCorr = max(max(corr));
       end
       
    end
    Coefs(bestWord,i) = 1;
    
end