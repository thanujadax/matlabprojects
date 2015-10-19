function sumstd = getSumStd(stdMat,dim)

varMat = stdMat.^2;

sumstd = sqrt(mean(varMat,dim));