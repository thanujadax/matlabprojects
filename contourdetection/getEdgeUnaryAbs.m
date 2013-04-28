function edgeUnary = getEdgeUnaryAbs(edgepixels,absOFR)
% unary value recalculation

% for each edge, get the absolute max OFR value of its each pixel and take the average
% as a unary score to be assigned to that edge

[numEdges,~] = size(edgepixels);
edgeUnary = zeros(numEdges,1);

for i=1:numEdges
	edgepix_i = edgepixels(i,:);
	edgepix_i = edgepix_i(edgepix_i>0);
	meanAbsOFR_i = mean(absOFR(edgepix_i));
	edgeUnary(i) = meanAbsOFR_i;		

end
