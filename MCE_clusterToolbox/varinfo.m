function dvi = varinfo(L1,L2)

[L1m,L1n] = size(L1);
[L2m,L2n] = size(L2);

numL1 = L1m*L1n;
numL2 = L2m*L2n;

L1 = reshape(L1,numL1,1);
L2 = reshape(L2,numL2,1);

%VARINFO Variation of information.
% DVI = VARINFO(L1,L2) returns the variation of information shared by two
% N-by-1 integer arrays of classification data, L1 and L2.
%
% Copyright (2009) Sandia Corporation. Under the terms of Contract 
% DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
% certain rights in this software.

% entropies for individual classification arrays
h1 = infoentropy(L1);
h2 = infoentropy(L2);

% mutual information shared by the classification arrays
i12 = mutualinfo(L1,L2);

% variation of inofrmation
dvi = h1 + h2 - 2*i12;