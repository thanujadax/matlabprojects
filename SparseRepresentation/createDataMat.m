newImagePath = '/home/thanuja/Dropbox/data/em_2013january/membranes/01.tif';
bb = 16;
slidingDist = 1;

IMin2 = imread(newImagePath);
IMin2 = IMin2(1:512,1:512);
[numRows numCols] = size(IMin2);
% normalize
IMin2 = IMin2./(max(max(IMin2)));

% using bb and slidingDist, generate label values compatible with the
% sparse codes
% i.e for bb=16, take the mode label of the 4 center pixels

numBlocks = prod(floor((size(IMin2)-bb)/slidingDist)+1);

labelVec = zeros(numBlocks,1);

startRow = bb/2;
endRow = numRows - bb/2;
startCol = bb/2;
endCol = numCols - bb/2;

j=0;
for r=startRow:slidingDist:endRow
    for c=startCol:slidingDist:endCol
        j=j+1;
        % get the mode label of the 4 pixels around j
        midMat = IMin2(r:r+1,c:c+1);
        membranes = find(midMat);
        if (membranes>2)
            labelVec(j)=1;
        end
    end
end
