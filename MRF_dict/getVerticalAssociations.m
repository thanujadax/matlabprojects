function AMat = getVerticalAssociations(coefmat,rowSize)

% coefmat - sparse matrix of coefficients. Rows correspond to the
% dictionary words and the columns correspond to the overlapping patches of
% the image.
% rowSize - number of patches in one row of the image


AMat = ones(size(coefmat,1));     % initializing association matrix
    
for i = 1:(size(coefmat,2)-rowSize)
    ind1 = find(coefmat(:,i));
    j = i + rowSize;            % downward neighbor of patch i
    ind2 = find(coefmat(:,j));  % index of the word for the next patch
    if (size(ind2,1)>0 && size(ind1,1)>0)
        AMat(ind1,ind2) = AMat(ind1,ind2) + 1;
    end
end

% save matrix
% save(outfilename,AMat);