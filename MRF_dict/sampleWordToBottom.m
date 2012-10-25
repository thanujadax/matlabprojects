function [wordIndV,probV] = sampleWordToBottom(wordToTop,Vmat)

% Inputs:
% Vmat - vertical associations given in element (top,bottom)
% wordToLeft - word to left, to be conditioned up on

% Outputs:
% wordIndV - sampled word
% probV - with probability

%%
% obtain associations for the relevant word
wordsToBottom = Hmat(wordToTop,:);
sumWordsToBottom = sum(wordsToBottom);
wordsToBottomNormalized = wordsToBottom./sumWordsToBottom; % contains the normalized cond. prob.

wordIndV = 0;
probV = 0;
% sampling
cumulativeProb = 0;          % stores the cumulative probability
x = rand(1);
for i = 1:length(wordsToBottomNormalized)
    cumulativeProb = cumulativeProb + wordsToBottomNormalized(i);
    if (cumulativeProb > x)
        wordIndV = i;
        probV = wordsToBottomNormalized(i);
        break;
    end
end

if (wordIndV == 0)
    wordIndV = length(wordsToBottomNormalized);
    probV = wordsToBottomNormalized(wordIndV);
end

