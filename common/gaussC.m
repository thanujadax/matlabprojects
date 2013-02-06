function val = gaussC(y, x, sigma, center)
xc = center(1);
yc = center(2);
% exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
exponent = ((x-xc).^2)./(2*sigma(1)) + ((y-yc).^2)./(2*sigma(2));
val       = (exp(-exponent));