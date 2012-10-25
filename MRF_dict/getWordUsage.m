% function [cumulative, usage] = getWordUsage(sparsecoefmatrix)

% generate statistics for each word based on their usage for image
% reconstruction

% outputs:
% usageSum = cumulative value of all the (positive) coefficients for
% each word. column vector.
% usageNum = number of times each word has been used. column vector.

function [usageNumMat, usageSumMat] = getWordUsage(sparsecoefmatrix)
cumulative = sum(sparsecoefmatrix,1); % sums up the coefficient values for each word

usageNumMat = zeros(size(sparsecoefmatrix,1),1);
for i = 1:size(sparsecoef,1)
    usageNumMat = length(find(sparsecoef(i,:))); % 
end

figure(1);
hist(cumulative,length(cumulative));
title('Cumulative coefficient values for each word');

figure(2);
hist(usage,length(usage));
title('Usage for each word');