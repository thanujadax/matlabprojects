function x0 = getInitValues(numEdges,numJ3,numJ4)
numVariables = numEdges*2 + numJ3*4 + numJ4*7;
x0 = ones(numVariables,1);
% edges - set the off states to zero
i=1:2:(numEdges*2-1);
x0(i) = 0;
% J3 set the off states to zero
j3Start = numEdges*2+1;
i = j3Start:4:(j3Start+4*numJ3-1);
x0(i) = 0;
% J4 set the off states to zero
j4Start = numEdges*2 + numJ3*4 + 1;
i = j4Start:7:(j4Start+7*numJ4-1);
x0(i) = 0;