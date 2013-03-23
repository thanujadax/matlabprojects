function dTheta = getAngleDifferences(angles0)
% for the angles given in columns 1, 2, 3 etc, return the differences in
% the order (theta1-theta2), (theta1-theta3), (theta2-theta3) in each
% column
% for nodes with any degree (i.e. for any number of columns in 'angles' the
% necessary combinations are calculated

if(numel(angles0)==1 && angles0(1)==0)
    dTheta = -1;
else
    % angles = angles0(angles0>0);
    [r,c] = find(angles0>0);
    angles(r,c) = angles0(r,c);
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
        dTheta = -1;  % no angles for this junction type
    end
end
