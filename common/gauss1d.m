function y = gauss1d(x,mean,sigma)

exponent = -1 .* ((x-mean).^2)./(2*(sigma.^2));

% y = (1/(sigma*(sqrt(2*pi)))).*exp(exponent);
y = exp(exponent);