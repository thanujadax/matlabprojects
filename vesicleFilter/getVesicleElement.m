function vesicleElement3D = getVesicleElement(innerRadius,outerRadius)
% returns the vesicle element to be used in the 3D convolution to detect
% vesicle like structures

vesicleElement3D = zeros(outerRadius*2,outerRadius*2,outerRadius*2);
[x y z] = meshgrid(1:outerRadius*2,1:outerRadius*2,1:outerRadius*2);

% create outerSphere
