function f1 = getUnaryPotentialVec(data,Dictionary)

numWords = size(Dictionary,2);

f1 = zeros(numWords,1);         % initialize column vector

dataAsColVector = repmat(data,numWords,1);

f1 = sqrt(sum((Dictionary' - dataAsColVector).^2),2);