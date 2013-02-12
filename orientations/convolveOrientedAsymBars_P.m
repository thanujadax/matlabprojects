function orientationScoreSpace3D = convolveOrientedAsymBars_P(img,barLength,barWidth,orientations,negLines)

% inputs: 
%   img - contains the pixels that should be used for calculating the hough
%   space
%   thresh - only the pixels above this values would be taken into account

% output:
% orientationScoreSpace3D

% each point in concern, votes for all the bars that it is contributing, at once
[numRows numCols] = size(img);
numOrientations = size(orientations,2);

% create the kernel
bar = ones(barWidth,barLength);

% twin shaded kernel
% halfW = floor(barWidth/2);
% bar(1:halfW,:) = -1;
negLine = ones(negLines,barLength).*-1;
bar = [bar;negLine];

% init
orientationScoreSpace3D = zeros(numRows,numCols,numOrientations);

parfor i=1:numOrientations
    img_i = img;
    % get the gaussian kernel
    bar_i = bar;
    % rotate appropriately
    orientation = orientations(i);
    if(orientation>0)
        bar_i = imrotate(bar_i,orientation);
    end
    % convolve to get the orientation score map for this orientation i
    % convResult = conv2(img_i,bar_i);
    convResult = xcorr2(img_i,bar_i);
    [convRows convCols] = size(convResult);
    rowMargin = (convRows-numRows)/2;
    colMargin = (convCols-numCols)/2;
    startR = rowMargin+1;
    stopR = convRows-rowMargin;
    startC = colMargin+1;
    stopC = convCols-colMargin;
    orientationScoreSpace3D(:,:,i) = convResult(startR:stopR,startC:stopC);
  % progressbar(j/totPoints); % update progress bar
end
% remove the border of size barWidth from the score
margin = barWidth+1;
orientationScoreSpace3D(1:margin,:,:) = 0;
orientationScoreSpace3D((numRows-margin):numRows,:,:) = 0;
orientationScoreSpace3D(:,1:margin,:) = 0;
orientationScoreSpace3D(:,(numCols-margin):numCols,:) = 0;

% normalize the score
orientationScoreSpace3D = orientationScoreSpace3D./...
        (max(max(max(orientationScoreSpace3D))));