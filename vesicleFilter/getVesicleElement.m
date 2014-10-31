function vesicleElement3D = getVesicleElement(innerRadius,outerRadius)
% returns the vesicle element to be used in the 3D convolution to detect
% vesicle like structures

vesicleElement3D = zeros(outerRadius*2,outerRadius*2,outerRadius*2);
[x y z] = meshgrid(1:outerRadius*2,1:outerRadius*2,1:outerRadius*2);

% center coordinates
px = outerRadius;
py = outerRadius;
pz = outerRadius;

% create outerSphere
inds_outerSphere = find((x-px).^2 + (y-py).^2 + (z-pz).^2 <= outerRadius.^2);
inds_innerSphere = find((x-px).^2 + (y-py).^2 + (z-pz).^2 <= innerRadius.^2);

vesicleElement3D(inds_outerSphere) = -1;
vesicleElement3D(inds_innerSphere) = 1;



