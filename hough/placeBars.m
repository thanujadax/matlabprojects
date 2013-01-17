function output = placeBars(output,goodpoints,peaks3D,orientations,barLength,barWidth)
% places bars at the given points at the right orientations and assigns
% them the right values

% Inputs:
%   goodpoints - indices of the points that have to be plotted with bars on
%   top of them
%   peaks3D - 3D hough space

[numRows numCols numOrientations] = size(peaks3D);
% goodpoints have 3 coordinates [row col ori]
[r,c,ori] = ind2sub([numRows,numCols,numOrientations],goodpoints);

vote = peaks3D(goodpoints);
output_vote = output(:,:,3); % value
output_ori = output(:,:,1); % hue    
progressbar('placing bars...')
numpoints = numel(goodpoints);
for i = 1:numpoints
   barInd = getBarPixInd(r(i),c(i),orientations(ori(i)),...
            barLength,barWidth,numRows,numCols);
   % update only the pixels that have a higher vote than the current vote
   barPixIndToUpdate = find(output_vote(barInd)<vote(i));
   output_vote(barInd(barPixIndToUpdate)) = vote(i);
   output_ori(barInd(barPixIndToUpdate)) = orientations(ori(i))/180;
   progressbar(i/numpoints);
end

output(:,:,1) = output_ori;
output(:,:,3) = output_vote;