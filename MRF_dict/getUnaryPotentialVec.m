function f1 = getUnaryPotentialVec(data,Dictionary,sigma)
% data is a col vector

numWords = size(Dictionary,2);

% f1 = zeros(numWords,1);         % initialize column vector

dataAsColVector = repmat(data',numWords,1);

Dictionary = Dictionary - min(min(Dictionary));
Dictionary = Dictionary ./ max(max(Dictionary));
Dictionary = Dictionary.*255;

% f1 = sqrt(sum(((Dictionary' - dataAsColVector).^2),2)/(2*sigma^2 * size(Dictionary,1)));

f1 = sqrt(sum(((Dictionary' - dataAsColVector).^2),2)/(size(Dictionary,1)));