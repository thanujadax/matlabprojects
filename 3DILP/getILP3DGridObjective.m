function f = getILP3DGridObjective(W,unaryScoresMat)

numVar = numel(unaryScoresMat);
[numRows,numCols] = size(unaryScoresMat);

f = zeros(numVar,1);

stop = 0;
for i=1:numCols
    start = stop +1;
    stop = stop + numRows;
    f(start:stop) = unaryScoresMat(:,i) .* W(i);
end
