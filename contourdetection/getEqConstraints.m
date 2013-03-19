function [Aeq,beq] = getEqConstraints(numEdges,j3Edges,j4Edges)
numJ3 = size(j3Edges,1);
numJ4 = size(j4Edges,1);
numCols = numEdges*2 + numJ3*4 + numJ4*7;
numRows = numEdges + numJ3*2 + numJ4*2;

Aeq = zeros(numRows,numCols);
beq = zeros(numRows,1);

beq(1:(numEdges+numJ3+numJ4)) = 1;

%% activation/inactivation constraints for each edge
j = 1;
for i=1:numEdges
    Aeq(i,j:(j+1)) = 1;
    j = j+2;
end

%% activation/inactivation constraints for each junction node
% J3
j = numEdges*2+1;
for i=(numEdges+1):(numEdges+numJ3)
    Aeq(i,j:(j+3)) = 1;
    j = j+4;
end
% J4
j = numEdges*2+numJ3*4+1;
for i=(numEdges+numJ3+1):(numEdges+numJ3+numJ4)
    Aeq(i,j:(j+6)) = 1;
    j = j+7;
end

%% closedness constraints
% for all nodes, get the active edge state variable indices
j3ActiveEdgeColInds = j3Edges .*2;
j4ActiveEdgeColInds = j4Edges .*2;

% J3
jColId = numEdges*2;
k = 1;
for i=(numEdges+numJ3+numJ4+1):(numEdges+numJ3*2+numJ4)
   Aeq(i,j3ActiveEdgeColInds(k,:)) = 1;          % marks active edges for junction i
   k = k+1;
   jIds = (jColId+2):(jColId+4);        % activate the 3 active states coeff for J3
   Aeq(i,jIds) = -2;                      % refer closedness constraint formulation
   
   % finally
   jColId = jColId + 4;
end

% J4
jColId = numEdges*2 + numJ3*4;
k = 1;
for i=(numEdges+numJ3*2+numJ4+1):(numEdges+numJ3*2+numJ4*2)
   Aeq(i,j4ActiveEdgeColInds(k,:)) = 1;          % marks active edges for junction i
   k = k+1;
   jIds = (jColId+2):(jColId+7);        % activate the 6 active states coeff for J4
   Aeq(i,jIds) = -2;                      % refer closedness constraint formulation
   
   % finally
   jColId = jColId + 7;
end

