function gradMap = getGradientMap(imgInv,sigma)
% returns a two layered matrix with the first layer being the grad
% magnitude and the second layer being the orientations scaled from 0 to
% 180

[dx, dy] = smoothGradient(imgInv, sigma);
% Calculate Magnitude of Gradient
magGrad = hypot(dx, dy);
% Normalize
magmax = max(magGrad(:));
if magmax > 0
    magGrad = magGrad / magmax;
end
angles = atan2(dy,dx).*180./pi;
angles = (angles + 180)./2;  % scaled from 0 to 180

gradMap(:,:,1) = magGrad;
gradMap(:,:,2) = angles;