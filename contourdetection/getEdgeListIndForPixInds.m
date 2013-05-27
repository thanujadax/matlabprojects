function edgeListInd = getEdgeListIndForPixInds(pixInds,edgePixels)
% Inputs
%   pixInds - the pix indices for which we want to find the edgeListInds
%   edgePixels - rowID is the edgeListInd. Each row contains the pixInds
%   for the corresponding edge

numPix = numel(pixInds);
k=0;    % iterator for the identified edges list
for i=1:numPix
    % get the corresponding edgeListInd for the pixels
    clear edgeListInd_i;
    [edgeListInd_i,~] = find(edgePixels==pixInds(i));
    if(~isempty(edgeListInd_i))
        k = k+1;
        edgeListInd(k) = edgeListInd_i;       
    end
end
edgeListInd = unique(edgeListInd);