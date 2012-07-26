function outMat = simpleThreshold(inMat, threshold)
% Simple thresholding function
% Every pixel < threshold is marked 0 and
% every pixel >= threshold is marked 1

sizeInMat = size(inMat);
outMat = zeros(sizeInMat);
for i=1:sizeInMat(1)
    for j=1:sizeInMat(2)
        if (inMat(i,j) >= threshold)
            outMat(i,j) = 1;
        end
    end 
end