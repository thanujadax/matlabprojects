function AMat = getHorizontalAssociations(coefmat,rowSize)

% coefmat - sparse matrix of coefficients. Rows correspond to the
% dictionary words and the columns correspond to the overlapping patches of
% the image.

AMat = ones(size(coefmat,1));       % initializing association matrix

ind1 = find(coefmat(:,1));          % index of the nonzero element of column 1 (index of atom for patch #1)
    
for i = 2:(size(coefmat,2))
    ind2 = find(coefmat(:,i));      % index of the word for the next patch
    if(mod(i,rowSize)==1)
        % the patch i is the first word of the row. not to the right of the
        % one before. do not count.
    else
    AMat(ind1,ind2) = AMat(ind1,ind2) + 1;  
    end
    ind1 = ind2;
end


