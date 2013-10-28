function [meanA,maxA,minA,stdMean,stdMax,stdMin,...
        medianMean,medianMax,medianMin] = getMatStats(A)

meanA = mean(mean(A));
maxA = max(max(A));
minA = min(min(A));

sd = std(A); % std vector containing std for each dimension
stdMax = max(sd);
stdMin = min(sd);
stdMean = mean(sd);

md = median(A); % median of each dimension
medianMax = max(md);
medianMin = min(md);
medianMean = mean(md);