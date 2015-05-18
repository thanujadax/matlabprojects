function edges2pixels = getEdges2Pixels(edgePixLabels)
% returns a N-by-(n+1) array of all the edges and their corresponding pixels.
% edgeID is given by the first column

% Input
%   edgePixLabels - N-by-2 array. each row is for one pixel with its edge
%   label in the 2nd col

%   pEdges2nodes - each row has two nodeInds

% numPixels = size(edgePixLabels,2);
numEdges = max(edgePixLabels(:,2));
edges2pixels = zeros(numEdges,1);

for i=1:numEdges
    pixInd = find(edgePixLabels(:,2)==i);
    if(~isempty(pixInd))
        pixInd = edgePixLabels(pixInd,1);
        % edges2pixels(i) = pixInd;
        if(~isempty(pixInd))
            edges2pixels(i,1)=i;        % looks trivial but will be useful when self edges are removed
            for j=1:numel(pixInd)
                edges2pixels(i,(j+1))=pixInd(j);
            end
        end
    end
end