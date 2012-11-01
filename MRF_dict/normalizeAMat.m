% to normalize the association matrices which contain counts
function normalizedMat = normalizeAMat(inMat,normalizeColumns)

% inputs
% inMat: matrix to be normalized
% normalizeColumns: if set, normalization is done columnwise. if zero,
% normalization is done per each row.

normalizedMat = zeros(size(inMat));

if(normalizeColumns)
    colSum = sum(inMat,1);      % gets a row vector containing the columnwise sums

    for i = 1 : size(inMat,2)
        normalizedMat(:,i) = inMat(:,i)./colSum(i);      % column i normalized
    end
else
        rowSum = sum(inMat,2);      % gets a col vector containing the columnwise sums

    for i = 1 : size(inMat,1)
        normalizedMat(i,:) = inMat(i,:)./rowSum(i);      % row i normalized
    end

end
