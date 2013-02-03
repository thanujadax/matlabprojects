function [output RGBimg] = reconstructHSVbars(peaks3D,orientations,threshFrac)
% output is a 3 layered matrix [H S V]
%   H - orientation
%   S - number of equally high orientations
%   V - vote
[numRows numCols numOrientations] = size(peaks3D);
maxVote = max(max(max(peaks3D)));
thresh = maxVote*threshFrac;
threshVotes3D = zeros(numRows,numCols,numOrientations);
indVote = peaks3D>thresh;
threshVotes3D(indVote) = peaks3D(indVote);
% construct the final output so that each pixel displays its maximum vote
display('Assembling the final output...');
progressbar('Assembling the final output')
output = zeros(numRows,numCols,3);
for r=1:numRows
    for c=1:numCols
        vote=max(threshVotes3D(r,c,:)); % vote stored in layer 3 - H
        if(vote==0)
            continue;
        else
            %output(r,c,3) = vote/maxVote;  % stores normalized vote
            output(r,c,3) = 0.5; % V - stores only if the point should be displayed
            output(r,c,2) = 1; % S
            ori = find(threshVotes3D(r,c,:)==vote);
            if(numel(ori)>1)
                ori = ori(1);
                output(r,c,2)=numel(ori)/numOrientations;
            elseif(numel(ori)==1)
                ori = orientations(ori);
            end
            output(r,c,1) = ori/180;           % H - orientation
        end
    end
    progressbar(r/numRows);
end

% create HSV image
hsvImage = cat(3,output(:,:,1),output(:,:,2),output(:,:,3));
% convert it to an RGB image
RGBimg = hsv2rgb(hsvImage);
