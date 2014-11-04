function vesicleElement2D = getVesicleElement2D(r1,r2,r3)

% r1 < r2 < r3

vesicleElement2D = zeros(r3*2,r3*2);
[x y] = meshgrid(1:r3*2,1:r3*2);
% center coordinates
px = r3;
py = r3;

% create sphere r3 (outer most)
inds_sphere_r3 = find((x-px).^2 + (y-py).^2 <= r3.^2);
% create sphere r2 (middle)
inds_sphere_r2 = find((x-px).^2 + (y-py).^2 <= r2.^2);
% create sphere r1 (inner most)
inds_sphere_r1 = find((x-px).^2 + (y-py).^2 <= r1.^2);

% the middle sphere is dark
vesicleElement2D(inds_sphere_r3) = 1; % bright
vesicleElement2D(inds_sphere_r2) = -1; % dark
vesicleElement2D(inds_sphere_r1) = 1; % bright