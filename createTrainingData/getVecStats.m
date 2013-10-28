function [meanA,maxA,minA,sd,md] = getVecStats(A)

meanA = mean(A);
maxA = max(A);
minA = min(A);

sd = std(A); % std vector containing std for each dimension
md = median(A); % median of each dimension