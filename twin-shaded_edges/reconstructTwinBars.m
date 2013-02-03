function output = reconstructTwinBars(peaks3D,orientations,barLength,barWidth)
% reconstructs the original image based on the peaks given for each
% orientation using the twin shaded bars

% Parallel

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

% output is a 3 layered matrix [H S V]
%   H - orientation
%   S - number of equally high orientations
%   V - vote

[numRows numCols numOrientations] = size(peaks3D);
outputs = zeros(numRows,numCols,numOrientations);
%numBarPix = barLength*barWidth;
maxVote = max(max(max(peaks3D)));
thresh = maxVote*0.7;
display('calculating pixel value per each orientation');

parfor i=1:numOrientations
    ori = orientations(i);
    voteMat = peaks3D(:,:,i);
    peaksInd = find(voteMat>thresh);
    numPeaks = numel(peaksInd);
    %barInd = zeros(numPeaks,numBarPix);
    %barVote = zeros(numPeaks,1);
    output_tmp = zeros(numRows,numCols);
    for j=1:numPeaks
       % place a bar on each peak as described above
       % barInd = getBar(numRows,numCols,peaksInd(j),barLength,barWidth,orientation);
       [r,c] = ind2sub([numRows,numCols],peaksInd(j));
       barInd = getBarPixInd(r,c,ori,barLength,barWidth,numRows,numCols);

       % a row of barInd corresponds to the peak j of this orientation   
       barVote = voteMat(peaksInd(j));
       % now barVote is sorted in the ascending order and barInd is
       % re-arranged to keep its correspondence with barVote

       % assign barVote to the barPixels only if the current value of
       % each barPixel is lower than barVote 
       % output_tmp(barInd) = barVote;   
       barPixIndToUpdate = find(output_tmp(barInd)<barVote);
       output_tmp(barInd(barPixIndToUpdate)) = barVote;
    end
    %[barVote barInd] = sortrows([barVote barInd],1);
    outputs(:,:,i) = output_tmp;    
end

% construct the final output so that each pixel displays its maximum vote
display('Assembling the final output...');
progressbar('Assembling the final output')
output = zeros(numRows,numCols,3);
for r=1:numRows
    for c=1:numCols
        vote=max(outputs(r,c,:)); % vote stored in layer 3 - H
        if(vote~=0)
            output(r,c,3) = vote/maxVote;
            ori = find(outputs(r,c,:)==vote);
            if(numel(ori)>1)
                ori = 0;
                output(r,c,2)=numel(ori)/numOrientations;
            else
                ori = orientations(ori);
            end
            output(r,c,1) = ori/180;           % orientation
        end
    end
    progressbar(r/numRows);
end
