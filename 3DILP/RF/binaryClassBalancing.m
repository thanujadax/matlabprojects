function [x,y] = binaryClassBalancing(fm,labels)

positiveLabels_logicalInds = (labels==1);
numPositiveLabels = sum(positiveLabels_logicalInds);

negativeLabels_logicalInds = (labels==0);
numNegativeLabels = sum(negativeLabels_logicalInds);

numBalancedDataPoints = min(numPositiveLabels,numNegativeLabels);

y_pos = labels(positiveLabels_logicalInds);
y_neg = labels(negativeLabels_logicalInds);

y = double([y_pos(1:numBalancedDataPoints); y_neg(1:numBalancedDataPoints)]);

x_pos = fm(positiveLabels_logicalInds,:);
x_neg = fm(negativeLabels_logicalInds,:);

x = [x_pos(1:numBalancedDataPoints,:); x_neg(1:numBalancedDataPoints,:)];