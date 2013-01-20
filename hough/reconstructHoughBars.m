function output = reconstructHoughBars(peaks3D,orientations,barLength,barWidth)
% reconstructs the original image based on the peaks given for each
% orientation

% Inputs:
%   peaks3D - a 3D array [row col orientation] containing the votes
%   orientations - e.g. [0 45 90 135]. can be any value from 1 to 180.
%   barLength - 
%   barWidth -


% for each orientation
%   get the peaks
%   place a bar on the peak according to the  given orientation
%   the intensity(color) of the bar is proportional to the vote of that
%   peak
%   the pixel intensities don't add up. instead, the higher vote gets
%   preference to decide its intensity when two overlapping bars compete
%   for the same pixel


[numRows numCols numOrientations] = size(peaks3D);
output = zeros(numRows,numCols);
margin = (max(barLength,barWidth)+1)/2 + 1;

for i=1:numOrientations
   orientation = orientations(i);
   voteMat = peaks3D(:,:,i);
   % make sure the margin is set to zero
    voteMat(1:margin,:) = 0;
    voteMat((numRows-margin):numRows,:) = 0;
    voteMat(:,1:margin) = 0;
    voteMat(:,(numCols-margin):numCols) = 0;
   % peaksInd = find(voteMat(margin:(numRows-margin),margin:(numCols-margin)));
   peaksInd = find(voteMat);
   numPeaks = numel(peaksInd);
   
   for j=1:numPeaks
       % place a bar on each peak as described above
       % barInd = getBar(numRows,numCols,peaksInd(j),barLength,barWidth,orientation);
       [r,c] = ind2sub([numRows,numCols],peaksInd(j));
%        if(r<margin||c<margin||r>=(numRows-margin)||c>=(numCols-margin))
%            continue;
%        end
       barInd = getBarPixInd(r,c,orientation,barLength,barWidth,numRows,numCols);
       if(barInd==-1)
           continue;
       end
       barVote = voteMat(peaksInd(j));
       % assign barVote to the barPixels only if the current value of
       % each barPixel is lower than barVote 
       % output(barInd) = barVote;   
       barPixIndToUpdate = find(output(barInd)<barVote);
       output(barInd(barPixIndToUpdate)) = barVote;
   end
       
end