function dTheta = getAngleDifferences(angles0)
% for the angles given in columns 1, 2, 3 etc, return the differences in
% the order (theta1-theta2), (theta1-theta3), (theta2-theta3) in each
% column
% for nodes with any degree (i.e. for any number of columns in 'angles' the
% necessary combinations are calculated

angles = angles0(angles0>0);
if(~isempty(angles))
    [numNodes,numTheta] = size(angles);
    edgeOrder = 1:numTheta;
    numCombinations = nchoosek(numTheta,2);
    combinationsArray = nchoosek(edgeOrder,2);

    dTheta = zeros(numNodes,numCombinations);

    for i=1:numCombinations
       dTheta(:,i) = abs(angles(:,combinationsArray(i,1)) - angles(:,combinationsArray(i,2)));
    end
else
    angles = -1;  % no angles for this junction type
end

