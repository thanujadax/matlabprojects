function edges2pixels = getEdges2Pixels(edgePixLabels)
% returns a N-by-n array of all the edges and their corresponding pixels.
% the row index corresponds to the edgeID

% Input
%   edgePixLabels - N-by-2 array. each row is for one pixel with its edge
%   label in the 2nd col

numPixels = size(edgePixLabels,2);
numEdges = max(edgePixLabels(:,2));
edges2pixels = zeros(numEdges,1);

for i=1:numEdges
    pixInd = find(edgePixLabels(:,2)==i);
    if(~isempty(pixInd))
        pixInd = edgePixLabels(pixInd,1);
        edges2pixels(i) = pixInd;
    end
end


