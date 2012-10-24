function AMat = getHorizontalAssociations(coefmat,outfilename)

% coefmat - sparse matrix of coefficients. Rows correspond to the
% dictionary words and the columns correspond to the overlapping patches of
% the image.

AMat = ones(size(coefmat,1));     % initializing association matrix

ind1 = find(coefmat(:,1));       % index of the nonzero element of column 1 (index of atom for patch #1)
% if size(ind1,1) == 0
%     ind1 = 0;
% end
    
for i = 2:(size(coefmat,2)-1)
    ind2 = find(coefmat(:,i));  % index of the word for the next patch
    if (size(ind2,1)>0 && size(ind1,1)>0)
        AMat(ind1,ind2) = AMat(ind1,ind2) + 1;
    end
    ind1 = ind2;
end


