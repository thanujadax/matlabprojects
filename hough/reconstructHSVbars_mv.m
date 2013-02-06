function [output RGBimg] = reconstructHSVbars_mv(peaks3D,orientations,barLength,barWidth,threshFrac)
% output is a 3 layered matrix [H S V]
%   H - orientation
%   S - number of equally high orientations
%   V - vote

[numRows numCols numOrientations] = size(peaks3D);
maxVote = max(max(max(peaks3D)));
thresh = maxVote*threshFrac;
output = zeros(numRows,numCols,3);
threshVotes3D = zeros(numRows,numCols,numOrientations);
indVote = peaks3D>thresh;
threshVotes3D(indVote) = peaks3D(indVote);
    
% construct the final output so that each pixel displays its maximum vote
display('Assembling the final output...');
progressbar('Assembling the final output')

for r=1:numRows
    for c=1:numCols
        vote=max(threshVotes3D(r,c,:)); % vote stored in layer 3 - H
        if(vote==0)
            continue;
        else   
            output(r,c,2) = 1; % S. default for all points above the threshold
            % get the weighted orientation as the vote. 
            % weighted according to the vote fraction           
            votes = zeros(1,numOrientations);
            votes(1,:) = threshVotes3D(r,c,:);          
            ori = find(votes==vote); % gets the orientation ind for max vote
            if(numel(ori)>1)
                ori = ori(1);
                % output(r,c,2)=numel(ori)/numOrientations;
                %disp('multiple max votes found for point %d,%d',r,c);
            end
            ori = orientations(ori);                        
            output(r,c,1) = ori/180;           % H - orientation
            % output(r,c,3) = vote/maxVote;       % V - score
            output(r,c,3) = 1;       % V - score
        end
    end
    progressbar(r/numRows);
end

% create HSV image
hsvImage = cat(3,output(:,:,1),output(:,:,2),output(:,:,3));
% convert it to an RGB image
RGBimg = hsv2rgb(hsvImage);
figure;imshow(RGBimg)
