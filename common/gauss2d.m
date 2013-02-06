function mat = gauss2d(mat, sigma, center)
% inputs
%   mat - 
%   sigma - [sigmaX sigmaY]
%   center - [centerX centerY]
gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, sigma, center);