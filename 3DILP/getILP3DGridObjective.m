function f = getILP3DGridObjective(W,unaryScoresMat)

numVar = numel(unaryScoresMat);
[numRows,numCols] = size(unaryScoresMat);

f = zeros(numVar,1);

for i=1:numCols
   unaryScoresMat(:,i) = unaryScoresMat(:,i) .* W(i); 
end

stop = 0;
for i=1:numRows
    start = stop +1;
    stop = stop + numCols;
    f(start:stop) = unaryScoresMat(i,:);
end
