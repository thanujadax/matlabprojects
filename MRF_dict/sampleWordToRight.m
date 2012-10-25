function [wordIndH,probH] = sampleWordToRight(wordToLeft,Hmat)

% Inputs:
% Hmat - horizontal associations given in element (left,right)
% wordToLeft - word to left, to be conditioned up on

% Outputs:
% wordIndH - sampled word
% probH - with probability

%%
% obtain associations for the relevant word
wordsToRight = Hmat(wordToLeft,:);
sumWordsToRight = sum(wordsToRight);
wordsToRightNormalized = wordsToRight./sumWordsToRight; % contains the normalized cond. prob.

wordIndH = 0;
probH = 0;
% sampling
cumulativeProb = 0;          % stores the cumulative probability
x = rand(1);
for i = 1:length(wordsToRightNormalized)
    cumulativeProb = cumulativeProb + wordsToRightNormalized(i);
    if (cumulativeProb > x)
        wordIndH = i;
        probH = wordsToRightNormalized(i);
        break;
    end
end

if (wordIndH == 0)
    wordIndH = length(wordsToRightNormalized);
    probH = wordsToRightNormalized(wordIndH);
end



