function [output RGBimg] = reconstructHSVlines(peaks3D,orientations,maxLength,barWidth,threshFrac)
% output is a 3 layered matrix [H S V]
%   H - orientation
%   S - number of equally high orientations
%   V - vote
% places a line at each peak whose length is proportional to the vote and
% the color is determined by the orientation

[numRows numCols numOrientations] = size(peaks3D);
maxVote = max(max(max(peaks3D)));
thresh = maxVote*threshFrac;
output = zeros(numRows,numCols,3);

% barWidth = 1;

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
       % a row of barInd corresponds to the peak j of this orientation   
       barVote = voteMat(peaksInd(j));
       %barLength = ceil(maxLength*barVote/maxVote);
       barLength = maxLength;
       barInd = getBarPixInd(r,c,ori,barLength,barWidth,numRows,numCols);

       
       % now barVote is sorted in the ascending order and barInd is
       % re-arranged to keep its correspondence with barVote

       % assign barVote to the barPixels only if the current value of
       % each barPixel is lower than barVote 
       % output_tmp(barInd) = barVote;   
       barPixIndToUpdate = find(output_tmp(barInd)<barVote);
       output_tmp(barInd(barPixIndToUpdate)) = barVote;
    end
    %[barVote barInd] = sortrows([barVote barInd],1);
    threshVotes3D(:,:,i) = output_tmp;    
end
    
% construct the final output so that each pixel displays its maximum vote
display('Assembling the final output...');
progressbar('Assembling the final output')

for r=1:numRows
    for c=1:numCols
        vote=max(threshVotes3D(r,c,:)); % vote stored in layer 3 - H
        if(vote==0)
            continue;
        else
            %output(r,c,3) = vote/maxVote;  % stores normalized vote
            output(r,c,3) = 1; % V - stores only if the point should be displayed
            output(r,c,2) = 1; % S
            
            
%             ori = find(threshVotes3D(r,c,:)==vote);
%             if(numel(ori)>1)
%                 ori = ori(1);
%                 output(r,c,2)=numel(ori)/numOrientations;
%             elseif(numel(ori)==1)
%                 ori = orientations(ori);
%             end
%             output(r,c,1) = ori/180;           % H - orientation
        
            votes = zeros(1,numOrientations);
            votes(1,:) = threshVotes3D(r,c,:);
            totVote = sum(votes);
            w_votes = votes/totVote; 
            w_orientations = orientations/180;
            w_ori = w_votes.*w_orientations;
            w_ori = sum(w_ori);
            
            output(r,c,1) = w_ori;           % H - orientation
        
        end
    end
    progressbar(r/numRows);
end

% create HSV image
hsvImage = cat(3,output(:,:,1),output(:,:,2),output(:,:,3));
% convert it to an RGB image
RGBimg = hsv2rgb(hsvImage);
figure;imagesc(RGBimg)
