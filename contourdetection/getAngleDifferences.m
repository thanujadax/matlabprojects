function dTheta = getAngleDifferences(angles)
% for the angles given in columns 1, 2, 3 etc, return the differences in
% the order (theta1-theta2), (theta1-theta3), (theta2-theta3) in each
% column

[numNodes,numTheta] = size(angles);
numCols = nchoosek(numTheta,2);
dTheta = zeros(numNodes,numCols);
if (numTheta==3)
    % three angle case
    % t1 - t2
    dTheta(:,1) = abs(angles(:,1) - angles(:,2));
    % t1 - t3
    dTheta(:,2) = abs(angles(:,1) - angles(:,3));
    % t2 - t3
    dTheta(:,3) = abs(angles(:,2) - angles(:,3));
elseif(numTheta==4)
    % four angle case
    % t1 - t2
    dTheta(:,1) = abs(angles(:,1) - angles(:,2));
    % t1 - t3
    dTheta(:,2) = abs(angles(:,1) - angles(:,3));
    % t1 - t4
    dTheta(:,3) = abs(angles(:,1) - angles(:,4));
    % t2 - t3
    dTheta(:,4) = abs(angles(:,2) - angles(:,3));
    % t2 - t4
    dTheta(:,5) = abs(angles(:,2) - angles(:,4));
    % t3 - t4
    dTheta(:,6) = abs(angles(:,3) - angles(:,4));
else
    % not defined for any other number of angles. return error.
    dTheta = -1;
    disp('Error:getAngleDifferences: number of edges per node is not 3 or 4');
end
